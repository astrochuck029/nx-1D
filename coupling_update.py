#!/usr/bin/env python3
"""
TMG / Simcenter NX 11 ORTHO thermochemical coupling
====================================================
Material model: T700/M21 carbon/epoxy, Tranchard et al.

Workflow for every coupling interval:
    TMG run
      -> binary TEMPF
      -> copy + TMG ASCII conversion
      -> this script
      -> next XML
      -> TMG run
      -> repeat

IMPORTANT:
  * TEMPF is element-temperature based. Temperatures are mapped by the actual
    XML element ID; no nodal averaging is performed.
  * NX mesh/collector names are NOT hard-coded. The physical thermochemical
    domain is detected from the XML and existing thermal-load group references.
  * The SetList collector/group used by the existing Q_dec load is preferred.
    If Q_dec does not exist, the largest thermal-load group overlapping the
    3-D solid mesh is selected. If no usable group exists, a new collector
    named ThermoCouplingDomain is created automatically.
  * 3d_mesh(1), 3d_mesh(2), ... names and their element counts are never
    assumed.
  * Heated_Face / Rear_Face / Insulated are not used to define the volume
    domain unless the XML itself clearly points to them as the solid domain.
  * Generated state sets are stored in ElementList only. Existing SetList
    collectors are preserved so existing BCs do not lose their group names.
  * The PowerShell driver keeps the CLI:
        --tempf --xml --output --state --dt

The thermochemical equations below are the Tranchard T700/M21 model used for
this coupling, including the two competitive reactions and the temperature-
dependent virgin/degraded properties.
"""

import argparse
import copy
import os
import re
import tempfile
import xml.etree.ElementTree as ET

import numpy as np


# ============================================================
# Configuration — Tranchard T700/M21
# ============================================================

R_GAS = 8.31446261815324       # J/(mol K)
DEFAULT_DT = 0.5               # s
N_STATES = 20
STATE_DA = 0.05

RHO_V = 1575.0                 # kg/m3
RHO_E = 1165.0                 # kg/m3
DRHO_TOTAL = RHO_V - RHO_E

# Effective emissivity model used in the previous coupling formulation:
#     epsilon = (1 - alpha) * epsilon_v + alpha * epsilon_e
# The current NX radiation definition uses 0.85 as its existing effective
# emissivity.  We retain 0.90 as the virgin value and 0.85 as the degraded
# value; change these two constants here if your validated material data use
# different endpoints.
EPSILON_V = 0.90
EPSILON_E = 0.85

# Table 3: log10(A), Ea, n, Kcat, contribution/fraction
LOG10_A1 = 0.8585
LOG10_A2 = 8.7268
A1 = 10.0 ** LOG10_A1
A2 = 10.0 ** LOG10_A2
Ea1 = 58.7107e3                # J/mol
Ea2 = 146.6740e3               # J/mol
n1 = 1.1125
n2 = 2.0825
KCAT2 = 10.0 ** 0.8744

F1 = 0.5540
F2 = 0.9479

# Table 5: kJ/kg -> J/kg
Q_DEC1 = 259.54e3              # J/kg
Q_DEC2 = -152.22e3             # J/kg

T_TABLE_C = np.arange(20.0, 701.0, 10.0)

SOLID_ELEMENT_TYPES = {
    "HEXA8", "HEXA20",
    "TET4", "TET10",
    "PENTA6", "PENTA15",
}

STATE_PREFIX = "3d_mesh_state_"
MAT_PREFIX = "Alpha_State_"
QDEC_PREFIX = "Q_dec_E"
AUTO_GROUP_NAME = "ThermoCouplingDomain"


# ============================================================
# Thermophysical properties
# ============================================================

def Cp_virgin(T_C):
    return 2.8773 * T_C + 687.31


def Cp_degraded(T_C):
    return (
        5.1327e-7 * T_C**3
        - 2.0761e-3 * T_C**2
        + 2.599 * T_C
        + 662.53
    )


def k_virgin_inplane(T_C):
    return 7.4675e-3 * T_C + 2.7811


def k_virgin_thru(T_C):
    return 1.1113e-3 * T_C + 0.61391


def k_degraded_inplane(T_C):
    return (
        -3.5481e-6 * T_C**2
        + 6.4898e-3 * T_C
        + 2.3005
    )


def k_degraded_thru(T_C):
    return (
        7.4228e-10 * T_C**3
        - 4.1903e-7 * T_C**2
        + 2.3397e-4 * T_C
        + 7.7211e-2
    )


def mixture_properties(T_C, alpha):
    """Return rho, Cp, kX, kY, kZ for the selected decomposition state."""
    alpha = float(np.clip(alpha, 0.0, 1.0))
    T_C = np.asarray(T_C, dtype=float)

    rho = RHO_V - DRHO_TOTAL * alpha
    cp = (1.0 - alpha) * Cp_virgin(T_C) + alpha * Cp_degraded(T_C)
    kip = (1.0 - alpha) * k_virgin_inplane(T_C) + alpha * k_degraded_inplane(T_C)
    kth = (1.0 - alpha) * k_virgin_thru(T_C) + alpha * k_degraded_thru(T_C)

    # Symmetric laminate: X/Y are identical for this ORTHO model.
    return rho, cp, kip, kip.copy(), kth


# ============================================================
# Pyrolysis-gas sensible heat sink (Tranchard Part I Eq. 18)
# ============================================================

def Cp_gas(T_C):
    """Specific heat of released pyrolysis gases [J/(kg K)]."""
    T_C = np.clip(np.asarray(T_C, dtype=float), 0.0, 1500.0)
    return (
        3.5977e-7 * T_C**3
        - 9.2485e-4 * T_C**2
        + 1.0610 * T_C
        + 1256.6
    )


def gas_enthalpy_rise(T_C, T0_C=30.0):
    """Integrated gas sensible enthalpy rise from T0 to local T [J/kg]."""
    T_C = np.asarray(T_C, dtype=float)

    def H(T):
        return (
            (3.5977e-7 / 4.0) * T**4
            - (9.2485e-4 / 3.0) * T**3
            + (1.0610 / 2.0) * T**2
            + 1256.6 * T
        )

    return np.maximum(H(T_C) - H(T0_C), 0.0)


# ============================================================
# Tranchard two-step competitive kinetics
# ============================================================

def reaction_rates(T_K, a1, a2):
    """Return d(alpha1)/dt and d(alpha2)/dt [1/s]."""
    T_K = np.maximum(np.asarray(T_K, dtype=float), 1.0)
    a1 = np.asarray(a1, dtype=float)
    a2 = np.asarray(a2, dtype=float)

    remaining = np.clip(1.0 - a1 - a2, 0.0, 1.0)

    r1 = A1 * remaining**n1 * np.exp(-Ea1 / (R_GAS * T_K))
    r2 = (
        A2
        * remaining**n2
        * np.exp(-Ea2 / (R_GAS * T_K))
        * (1.0 + KCAT2 * np.clip(a2, 0.0, 1.0))
    )
    return r1, r2


def rk4_step(T_K, a1, a2, dt):
    """Advance both competitive reaction extents by one RK4 step."""
    if dt <= 0.0:
        raise ValueError("dt must be > 0")

    k1a, k1b = reaction_rates(T_K, a1, a2)
    k2a, k2b = reaction_rates(
        T_K, a1 + 0.5 * dt * k1a, a2 + 0.5 * dt * k1b
    )
    k3a, k3b = reaction_rates(
        T_K, a1 + 0.5 * dt * k2a, a2 + 0.5 * dt * k2b
    )
    k4a, k4b = reaction_rates(
        T_K, a1 + dt * k3a, a2 + dt * k3b
    )

    new_a1 = a1 + dt / 6.0 * (k1a + 2.0 * k2a + 2.0 * k3a + k4a)
    new_a2 = a2 + dt / 6.0 * (k1b + 2.0 * k2b + 2.0 * k3b + k4b)

    new_a1 = np.clip(new_a1, 0.0, 1.0)
    new_a2 = np.clip(new_a2, 0.0, 1.0)

    # Enforce the competitive constraint without imposing a false relation
    # such as alpha2 <= alpha1.
    total = new_a1 + new_a2
    over = total > 1.0
    if np.any(over):
        scale = np.ones_like(total)
        scale[over] = 1.0 / total[over]
        new_a1 *= scale
        new_a2 *= scale

    return new_a1, new_a2


def decomposition_degree(a1, a2):
    """alpha = (rho_v-rho)/(rho_v-rho_e)."""
    return np.clip(F1 * a1 + F2 * a2, 0.0, 1.0)


def density_from_reactions(a1, a2):
    """rho = rho_v - Delta_rho1*a1 - Delta_rho2*a2."""
    drho1 = F1 * DRHO_TOTAL
    drho2 = F2 * DRHO_TOTAL
    rho = RHO_V - drho1 * a1 - drho2 * a2
    return np.clip(rho, RHO_E, RHO_V)


def q_decomposition_from_average_rates(T_C, r1_avg, r2_avg):
    """
    Net volumetric decomposition heat generation [W/m3].

    The coupling interval is finite.  Therefore the most consistent heat
    source for TMG over that interval is based on the average reaction rates
    implied by the RK4 change in reaction extents:

        r_i,avg = (alpha_i[n+1] - alpha_i[n]) / dt

    The released-gas sensible enthalpy sink is then evaluated at the local
    TEMPF temperature.
    """
    drho1 = F1 * DRHO_TOTAL
    drho2 = F2 * DRHO_TOTAL

    q_chem = Q_DEC1 * drho1 * r1_avg + Q_DEC2 * drho2 * r2_avg
    drho_g = drho1 * r1_avg + drho2 * r2_avg
    dh_gas = gas_enthalpy_rise(T_C, T0_C=30.0)

    return q_chem - drho_g * dh_gas


# ============================================================
# XML helpers
# ============================================================

def fmt(value):
    return f" {float(value):.7E}"


MAX_TMG_UID = 99999
MATERIAL_UID_START = 100
QDEC_UID_START = 10000

def all_uids(root, exclude_generated=False):
    """Return integer UIDs present in the XML.

    TMG/NX 11 uses a five-digit UID field in several places.  The old
    implementation used the largest UID in the entire XML when generating
    new materials and Q_dec loads.  Because each solution contains thousands
    of Q_dec loads, that made the UID grow every coupling step and eventually
    overflow the TMG five-digit limit (the observed late-run crash).
    """
    used = set()
    for elem in root.iter():
        if exclude_generated:
            if is_generated_material(elem) if elem.tag == "Material" else False:
                continue
            if elem.tag == "ThermalLoad" and is_generated_qdec(elem):
                continue
        value = elem.get("uid")
        if value is None:
            continue
        try:
            uid = int(value)
        except (TypeError, ValueError):
            continue
        if 0 <= uid <= MAX_TMG_UID:
            used.add(uid)
    return used


def max_uid(root):
    """Largest currently-used UID, restricted to TMG's five-digit range."""
    return max(all_uids(root), default=0)


def allocate_uid_block(root, count, preferred_start):
    """Allocate a contiguous UID block without exceeding five digits."""
    if count <= 0:
        raise ValueError("count must be > 0")

    used = all_uids(root, exclude_generated=True)
    start = max(1, int(preferred_start))

    # First try the preferred deterministic range.
    for candidate in (start, max_uid(root) + 1):
        candidate = max(1, candidate)
        end = candidate + count - 1
        if end > MAX_TMG_UID:
            continue
        block = range(candidate, end + 1)
        if all(uid not in used for uid in block):
            return candidate

    # Last-resort search for any free contiguous block.
    run = 0
    run_start = None
    for uid in range(1, MAX_TMG_UID + 1):
        if uid in used:
            run = 0
            run_start = None
            continue
        if run == 0:
            run_start = uid
        run += 1
        if run >= count:
            return run_start

    raise RuntimeError(
        f"Could not allocate {count} unique TMG UIDs below {MAX_TMG_UID}."
    )


def get_element_list(root):
    node = root.find("ElementList")
    if node is None:
        raise RuntimeError("Could not find <ElementList> in XML.")
    return node


def get_set_list(root):
    node = root.find("SetList")
    if node is None:
        raise RuntimeError("Could not find <SetList> in XML.")
    return node


def parse_element_record(e):
    text = (e.text or "").strip()
    if not text:
        raise RuntimeError("Found an <E> entry with no element data.")
    values = text.split()
    try:
        eid = int(values[0])
    except ValueError as exc:
        raise RuntimeError(f"Invalid element ID in <E>: {text!r}") from exc
    if eid <= 0:
        raise RuntimeError(f"Invalid element ID {eid}.")
    return eid, values


def is_state_set(s):
    return s.get("uname", "").startswith(STATE_PREFIX)


def is_generated_material(m):
    return m.get("uname", "").startswith(MAT_PREFIX)


def is_generated_qdec(load):
    uname = load.get("uname", "")
    return uname == "Q_dec" or uname.startswith(QDEC_PREFIX)


def solid_element_sets(element_list):
    result = []
    for s in element_list.findall("Set"):
        if s.get("elementType", "") in SOLID_ELEMENT_TYPES and s.findall("E"):
            result.append(s)
    return result


def collect_from_element_sets(sets):
    data = {}
    types = {}
    names = {}

    for s in sets:
        etype = s.get("elementType", "")
        uname = s.get("uname", "")
        for e in s.findall("E"):
            eid, values = parse_element_record(e)
            if eid in data:
                raise RuntimeError(
                    f"Duplicate solid element ID {eid} found in XML."
                )
            data[eid] = values
            types[eid] = etype
            names[eid] = uname

    if not data:
        raise RuntimeError("No solid elements found in XML.")
    return data, types, names


def setlist_selection_ids(s):
    """Read NX SetList collector IDs from <Selection><_>...</_></Selection>."""
    ids = []
    for sel in s.findall("Selection"):
        for item in sel.findall("_"):
            text = (item.text or "").strip()
            if not text:
                continue
            try:
                ids.append(int(text))
            except ValueError:
                continue
        # Some exports use <el> rather than <_> inside a Selection.
        for item in sel.findall("el"):
            text = (item.text or "").strip()
            if not text:
                continue
            try:
                ids.append(int(text))
            except ValueError:
                continue
    return sorted(set(ids))


def thermal_load_group_names(root):
    """Return group names referenced by thermal loads, preserving load order."""
    names = []
    loads = root.find("./Loads/ThermalLoadList")
    if loads is None:
        return names
    for load in loads.findall("ThermalLoad"):
        for sel in load.findall("Selection"):
            name = sel.get("groupName")
            if name and name not in names:
                names.append(name)
    return names


def group_name_from_generated_qdec(root):
    """Prefer the existing Q_dec selection group when it is still valid."""
    loads = root.find("./Loads/ThermalLoadList")
    if loads is None:
        return None
    for load in loads.findall("ThermalLoad"):
        if load.get("uname", "") != "Q_dec":
            continue
        for sel in load.findall("Selection"):
            name = sel.get("groupName")
            if name:
                return name
    return None


def build_setlist_inventory(root):
    inventory = {}
    for s in get_set_list(root).findall("Set"):
        name = s.get("uname", "")
        if name:
            inventory[name] = set(setlist_selection_ids(s))
    return inventory


def detect_physical_domain(root, all_element_data):
    """
    Detect the physical thermochemical solid domain.

    Priority:
      1. Existing generated state sets (subsequent coupling cycle).
      2. A SetList collector referenced by Q_dec and overlapping solids.
      3. The SetList collector with the largest overlap with the solid mesh
         among all thermal-load group references.
      4. The SetList collector with the largest overlap with the solid mesh.
      5. All solid element sets.

    This is intentionally based on element IDs, not names like 3d_mesh(1).
    """
    element_list = get_element_list(root)
    solid_sets = solid_element_sets(element_list)
    if not solid_sets:
        raise RuntimeError("No 3-D solid element sets found in XML.")

    solid_ids = set(all_element_data)

    # On subsequent cycles the generated state sets are the authoritative
    # physical mesh. They contain the exact same original EIDs.
    state_sets = [s for s in solid_sets if is_state_set(s)]
    if state_sets:
        data, types, names = collect_from_element_sets(state_sets)
        return data, types, names, "generated state sets"

    set_inventory = build_setlist_inventory(root)
    load_groups = thermal_load_group_names(root)
    q_group = group_name_from_generated_qdec(root)

    # Rank candidate collectors.  A group covering all solid IDs wins.
    candidates = []
    for name, ids in set_inventory.items():
        overlap = ids & solid_ids
        if not overlap:
            continue
        is_q = (q_group == name)
        referenced = (name in load_groups)
        candidates.append(
            (
                len(overlap) == len(solid_ids),
                is_q,
                referenced,
                len(overlap),
                name,
                overlap,
            )
        )

    if candidates:
        candidates.sort(reverse=True, key=lambda x: (x[0], x[1], x[2], x[3], x[4]))
        _, is_q, referenced, count, name, overlap = candidates[0]

        selected_sets = [s for s in solid_sets if any(
            parse_element_record(e)[0] in overlap for e in s.findall("E")
        )]
        selected_data, selected_types, selected_names = collect_from_element_sets_filtered(
            selected_sets, overlap
        )
        if selected_data:
            reason = "SetList collector"
            if is_q:
                reason += " referenced by existing Q_dec"
            elif referenced:
                reason += " referenced by thermal load"
            print(f"Physical group : {name}")
            print(f"Physical group elements : {len(overlap)}")
            print(f"Physical group source : {reason}")
            return selected_data, selected_types, selected_names, f"collector '{name}'"

    # Last-resort safe fallback: all 3-D solid elements.
    data, types, names = collect_from_element_sets(solid_sets)
    return data, types, names, "all 3-D solid element sets"


def collect_from_element_sets_filtered(sets, allowed_ids):
    data = {}
    types = {}
    names = {}
    allowed_ids = set(allowed_ids)

    for s in sets:
        etype = s.get("elementType", "")
        uname = s.get("uname", "")
        for e in s.findall("E"):
            eid, values = parse_element_record(e)
            if eid not in allowed_ids:
                continue
            if eid in data:
                raise RuntimeError(
                    f"Duplicate physical solid element ID {eid} while using a collector."
                )
            data[eid] = values
            types[eid] = etype
            names[eid] = uname

    if set(data) != allowed_ids:
        missing = sorted(allowed_ids - set(data))
        preview = ", ".join(map(str, missing[:20]))
        if len(missing) > 20:
            preview += ", ..."
        raise RuntimeError(
            "Collector references solid element IDs that cannot be found in "
            f"ElementList: {preview}"
        )
    return data, types, names


def ensure_qdec_group(root, element_ids):
    """
    Detect the groupName used for Q_dec.

    If an existing Q_dec group is valid, reuse it. Otherwise choose the best
    existing SetList collector. If none is suitable, create a new collector.
    """
    element_ids = set(element_ids)
    set_list = get_set_list(root)
    inventory = build_setlist_inventory(root)

    existing_q = group_name_from_generated_qdec(root)
    if existing_q and existing_q in inventory and element_ids.issubset(inventory[existing_q]):
        print(f"Q_dec group  : {existing_q} (from existing Q_dec)")
        return existing_q

    # Prefer a collector that contains all physical IDs.
    full = [name for name, ids in inventory.items() if element_ids.issubset(ids)]
    if full:
        # Prefer one referenced by a thermal load, then shortest/smallest.
        refs = set(thermal_load_group_names(root))
        full.sort(key=lambda n: (n not in refs, len(inventory[n]), n))
        name = full[0]
        print(f"Q_dec group  : {name} (auto-detected collector)")
        return name

    # Otherwise choose the collector with the largest physical overlap.
    ranked = []
    refs = set(thermal_load_group_names(root))
    for name, ids in inventory.items():
        overlap = len(ids & element_ids)
        if overlap:
            ranked.append((overlap == len(element_ids), overlap, name in refs, name))
    if ranked:
        ranked.sort(reverse=True)
        name = ranked[0][3]
        print(f"Q_dec group  : {name} (largest physical overlap)")
        return name

    # No usable group: create a valid NX SetList collector automatically.
    used_uids = []
    for s in set_list.findall("Set"):
        try:
            used_uids.append(int(s.get("uid", "0")))
        except ValueError:
            pass
    uid = max(used_uids, default=0) + 1

    name = AUTO_GROUP_NAME
    existing_names = set(inventory)
    counter = 1
    while name in existing_names:
        counter += 1
        name = f"{AUTO_GROUP_NAME}_{counter}"

    s = ET.SubElement(set_list, "Set", {"uid": str(uid), "uname": name})
    sel = ET.SubElement(s, "Selection", {"step": "1", "class": "Elements"})
    for eid in sorted(element_ids):
        ET.SubElement(sel, "_").text = str(eid)

    print(f"Q_dec group  : {name} (created automatically)")
    return name


# ============================================================
# TEMPF reader
# ============================================================

_NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")


def parse_temperature_pairs(filename):
    """Read ASCII TEMPF and retain the final temperature for each element ID."""
    values = {}
    with open(filename, "r", errors="ignore") as f:
        for line in f:
            tokens = _NUM_RE.findall(line)
            if len(tokens) < 2:
                continue
            nums = []
            for token in tokens:
                try:
                    nums.append(float(token))
                except ValueError:
                    pass

            # TEMPF records are normally: element_id temperature.
            # Ignore the -99999 time-header ID.
            for i in range(len(nums) - 1):
                ident = nums[i]
                temp = nums[i + 1]
                if not np.isfinite(ident) or not np.isfinite(temp):
                    continue
                iid = int(round(ident))
                if iid < 1 or abs(ident - iid) > 1e-10:
                    continue
                if -100.0 <= temp <= 5000.0:
                    values[iid] = temp
                    break

    if not values:
        raise RuntimeError(f"No usable element-temperature records found in TEMPF: {filename}")
    return values


def map_temperatures_to_elements(temp_values, element_data):
    eids = sorted(element_data)
    missing = [eid for eid in eids if eid not in temp_values]
    if missing:
        preview = ", ".join(map(str, missing[:30]))
        if len(missing) > 30:
            preview += ", ..."
        raise RuntimeError(
            "TEMPF/XML ELEMENT-ID MISMATCH\n"
            f"  XML physical elements : {len(eids)}\n"
            f"  Required ID range     : {eids[0]} .. {eids[-1]}\n"
            f"  TEMPF matched         : {len(eids) - len(missing)}\n"
            f"  Missing IDs           : {preview}"
        )

    T_C = np.asarray([temp_values[eid] for eid in eids], dtype=float)
    print(
        f"TEMPF mapping : direct EID mapping ({len(eids)} physical elements)"
    )
    return T_C


# ============================================================
# Persistent decomposition state
# ============================================================

def load_state(filename, element_ids):
    n = len(element_ids)
    if not os.path.exists(filename):
        print(f"No previous decomposition state. Starting alpha=0 for {n} elements.")
        return np.zeros(n), np.zeros(n)

    with np.load(filename, allow_pickle=False) as data:
        required = {"element_ids", "a1", "a2"}
        if not required.issubset(data.files):
            raise RuntimeError(
                f"Incompatible state file: {filename}. Delete it once before a new run."
            )
        old_ids = np.asarray(data["element_ids"], dtype=np.int64)
        a1 = np.asarray(data["a1"], dtype=float)
        a2 = np.asarray(data["a2"], dtype=float)

    ids = np.asarray(element_ids, dtype=np.int64)
    if old_ids.shape != ids.shape or not np.array_equal(old_ids, ids):
        raise RuntimeError(
            "Saved decomposition state does not match the current physical "
            "element IDs. Use a new --state file when the mesh changes."
        )
    if a1.shape != ids.shape or a2.shape != ids.shape:
        raise RuntimeError("Saved alpha arrays do not match the current mesh.")

    a1 = np.clip(a1, 0.0, 1.0)
    a2 = np.clip(a2, 0.0, 1.0)
    total = a1 + a2
    over = total > 1.0
    if np.any(over):
        a1[over] /= total[over]
        a2[over] /= total[over]
    return a1, a2


def save_state_atomic(filename, element_ids, a1, a2):
    directory = os.path.dirname(os.path.abspath(filename)) or "."
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".decomp_", suffix=".npz", dir=directory)
    os.close(fd)
    try:
        np.savez(
            tmp,
            element_ids=np.asarray(element_ids, dtype=np.int64),
            a1=np.asarray(a1, dtype=float),
            a2=np.asarray(a2, dtype=float),
        )
        os.replace(tmp, filename)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


# ============================================================
# Effective emissivity / NX radiation update
# ============================================================

def effective_emissivity(alpha):
    """
    Linear virgin/degraded mixture law:

        epsilon = (1-alpha)*epsilon_v + alpha*epsilon_e

    Alpha is clipped only for the emissivity calculation so the returned
    value remains physically bounded between the two endpoint values.
    """
    a = np.clip(np.asarray(alpha, dtype=float), 0.0, 1.0)
    return (1.0 - a) * EPSILON_V + a * EPSILON_E


def radiation_face_element_ids(radiation):
    """Return solid element IDs referenced by an NX SERadiation Selection."""
    ids = []
    for sel in radiation.findall("Selection"):
        # NX face selections are commonly: <fa>element_id face_number</fa>
        for item in sel.findall("fa"):
            toks = (item.text or "").split()
            if not toks:
                continue
            try:
                ids.append(int(toks[0]))
            except ValueError:
                continue

        # Some exports may use element selections instead.
        for item in sel.findall("el"):
            toks = (item.text or "").split()
            if not toks:
                continue
            try:
                ids.append(int(toks[0]))
            except ValueError:
                continue

        # And some SetList-style selections use <_>.
        for item in sel.findall("_"):
            toks = (item.text or "").split()
            if not toks:
                continue
            try:
                ids.append(int(toks[0]))
            except ValueError:
                continue

    return sorted(set(ids))


def update_effective_emissivity(root, element_ids, alpha):
    """
    Update every NX Simple Radiation to Environment object for which the
    selected face elements can be matched to the physical solid EIDs.

    NX stores Effective Emissivity as one scalar on a SERadiation object,
    while alpha is available per solid element. Therefore the scalar used by
    that radiation object is the arithmetic mean of the effective emissivity
    of its selected physical surface elements.

    Face selections such as <fa>4008 2</fa> refer to a face of solid element
    4008, so the first number is the physical solid EID used here.
    """
    rad_list = root.find("./Constraints/SERadiationList")
    if rad_list is None:
        print("Radiation update  : no SERadiationList found; skipped")
        return []

    alpha_by_eid = {int(eid): float(a) for eid, a in zip(element_ids, alpha)}
    physical_ids = set(alpha_by_eid)
    updated = []

    for radiation in rad_list.findall("SERadiation"):
        face_ids = radiation_face_element_ids(radiation)
        matched = [eid for eid in face_ids if eid in physical_ids]
        if not matched:
            continue

        eps_values = effective_emissivity(
            np.asarray([alpha_by_eid[eid] for eid in matched], dtype=float)
        )
        eps = float(np.mean(eps_values))
        eps = float(np.clip(eps, 0.0, 1.0))

        prop = None
        for p in radiation.findall("Property"):
            if p.get("name", "") == "Effective Emissivity":
                prop = p
                break

        if prop is None:
            prop = ET.Element("Property", {"name": "Effective Emissivity"})
            # Insert after View Factor when possible, otherwise append.
            children = list(radiation)
            insert_at = len(children)
            for i, child in enumerate(children):
                if child.tag == "Property" and child.get("name", "") == "View Factor":
                    insert_at = i + 1
                    break
            radiation.insert(insert_at, prop)

        for child in list(prop):
            prop.remove(child)
        ET.SubElement(prop, "Value").text = fmt(eps)

        uname = radiation.get("uname", "")
        updated.append((uname, len(matched), eps, float(np.mean(
            [alpha_by_eid[eid] for eid in matched]
        ))))

    if updated:
        for uname, n, eps, a_mean in updated:
            print(
                f"Radiation         : {uname or '<unnamed>'} | "
                f"faces={n} | alpha_mean={a_mean:.6e} | "
                f"epsilon_eff={eps:.7f}"
            )
    else:
        print("Radiation update  : no SERadiation surface matched physical EIDs")

    return updated


# ============================================================
# Material XML generation
# ============================================================

def replace_property_with_table(material, name, values):
    old = None
    for p in material.findall("Property"):
        if p.get("name") == name:
            old = p
            break

    new = ET.Element(
        "Property", {"name": name, "type": "XYTable", "x": "Temperature"}
    )
    for T, value in zip(T_TABLE_C, values):
        ET.SubElement(new, "_").text = f"{T: .7E}    {float(value): .7E}"

    if old is None:
        material.append(new)
    else:
        idx = list(material).index(old)
        material.remove(old)
        material.insert(idx, new)


def replace_constant_property(material, name, value):
    old = None
    for p in material.findall("Property"):
        if p.get("name") == name:
            old = p
            break

    if old is None:
        old = ET.SubElement(material, "Property", {"name": name, "type": "Constant"})
    else:
        old.set("type", "Constant")
        old.attrib.pop("x", None)
        for child in list(old):
            old.remove(child)

    ET.SubElement(old, "Value").text = fmt(value)


def find_ortho_template(root):
    material_list = root.find("MaterialList")
    if material_list is None:
        raise RuntimeError("Could not find <MaterialList> in XML.")

    generated = [m for m in material_list.findall("Material") if is_generated_material(m)]
    candidates = [
        m for m in material_list.findall("Material")
        if m.get("type", "").upper() == "ORTHO" and not is_generated_material(m)
    ]
    if not candidates:
        raise RuntimeError("No original ORTHO material found in XML.")

    # Prefer the material used by an original physical element set.
    element_list = get_element_list(root)
    material_uids = set()
    for s in solid_element_sets(element_list):
        uid = s.get("material")
        if uid:
            material_uids.add(uid)
    for m in candidates:
        if m.get("uid") in material_uids:
            return m
    return candidates[0]


def build_material_states(root):
    material_list = root.find("MaterialList")
    if material_list is None:
        raise RuntimeError("Could not find <MaterialList> in XML.")

    template = find_ortho_template(root)

    for m in list(material_list.findall("Material")):
        if is_generated_material(m):
            material_list.remove(m)

    # IMPORTANT: do NOT base material UIDs on Q_dec load UIDs.
    # Q_dec has thousands of generated loads and using their maximum UID
    # caused UID growth from solution to solution and the late-run TMG crash.
    base_uid = allocate_uid_block(root, N_STATES, MATERIAL_UID_START)
    state_uids = {}

    alpha_values = np.arange(1, N_STATES + 1, dtype=float) * STATE_DA
    for state, alpha in enumerate(alpha_values, start=1):
        rho, cp, kx, ky, kz = mixture_properties(T_TABLE_C, alpha)
        mat = copy.deepcopy(template)
        mat.set("uid", str(base_uid + state - 1))
        mat.set("uname", f"{MAT_PREFIX}{state:02d}")
        mat.set("type", "ORTHO")

        replace_constant_property(mat, "Mass Density", rho)
        replace_property_with_table(mat, "Specific Heat", cp)
        replace_property_with_table(mat, "Thermal Conductivity X", kx)
        replace_property_with_table(mat, "Thermal Conductivity Y", ky)
        replace_property_with_table(mat, "Thermal Conductivity Z", kz)

        material_list.append(mat)
        state_uids[state] = base_uid + state - 1

    return state_uids


# ============================================================
# State sets in ElementList
# ============================================================

def rebuild_state_sets(root, element_data, element_types, states, material_uids):
    element_list = get_element_list(root)

    # Remove every ElementList solid set. This removes both the original mesh
    # partitions and previously generated state sets. SetList collectors are
    # deliberately untouched.
    for s in list(element_list.findall("Set")):
        if s.get("elementType", "") in SOLID_ELEMENT_TYPES and s.findall("E"):
            element_list.remove(s)

    expected = sorted(element_data)

    for state in range(1, N_STATES + 1):
        by_type = {}
        for eid in expected:
            if int(states[eid]) == state:
                by_type.setdefault(element_types[eid], []).append(eid)

        for etype, eids in sorted(by_type.items()):
            suffix = "" if len(by_type) == 1 else f"_{etype}"
            s = ET.Element(
                "Set",
                {
                    "uname": f"{STATE_PREFIX}{state:02d}{suffix}",
                    "elementType": etype,
                    "material": str(material_uids[state]),
                    "opticalProperty": "0",
                },
            )
            for eid in eids:
                e = ET.SubElement(s, "E")
                e.text = " ".join(element_data[eid])
            element_list.append(s)

    # Verify one-and-only-one appearance of every physical solid EID.
    seen = []
    for s in element_list.findall("Set"):
        if s.get("elementType", "") not in SOLID_ELEMENT_TYPES:
            continue
        if not s.get("uname", "").startswith(STATE_PREFIX):
            continue
        for e in s.findall("E"):
            seen.append(parse_element_record(e)[0])

    if sorted(seen) != expected:
        raise RuntimeError(
            "State-set rebuild verification failed: physical solid EIDs changed."
        )


# ============================================================
# Q_dec generation
# ============================================================

def create_qdec_load(uid, element_id, q_value, group_name):
    load = ET.Element(
        "ThermalLoad",
        {
            "uid": str(uid),
            "uname": f"{QDEC_PREFIX}{element_id}",
            "type": "Heat Generation",
        },
    )
    ET.SubElement(load, "Description")

    p = ET.SubElement(load, "Property", {"name": "Heat Generation"})
    ET.SubElement(p, "Value").text = fmt(q_value)

    for name, value in [
        ("Control Heater", "0"),
        ("Thermostat", "-1"),
        ("Specify Layer to Apply to", "0"),
        ("Apply to", "0"),
        ("Layer Number", "1"),
    ]:
        p = ET.SubElement(load, "Property", {"name": name})
        ET.SubElement(p, "Value").text = value

    sel = ET.SubElement(
        load,
        "Selection",
        {"step": "1", "groupName": group_name},
    )
    ET.SubElement(sel, "el").text = str(element_id)
    return load


def update_qdec_loads(root, element_ids, qdec, group_name):
    loads = root.find("./Loads/ThermalLoadList")
    if loads is None:
        raise RuntimeError("Could not find <Loads>/<ThermalLoadList> in XML.")

    for load in list(loads.findall("ThermalLoad")):
        if is_generated_qdec(load):
            loads.remove(load)

    # Reuse a bounded UID block on every coupling step.  Old generated Q_dec
    # loads were removed above, so their previous UIDs must never drive the
    # next solution's UID upward.
    uid = allocate_uid_block(root, len(element_ids), QDEC_UID_START)
    for eid, q in zip(element_ids, qdec):
        loads.append(create_qdec_load(uid, eid, q, group_name))
        uid += 1


# ============================================================
# XML writer / verification
# ============================================================

def indent_xml(elem, level=0):
    i = "\n" + level * "  "
    j = "\n" + (level - 1) * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for child in elem:
            indent_xml(child, level + 1)
        if not elem[-1].tail or not elem[-1].tail.strip():
            elem[-1].tail = i
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = j


def write_xml_atomic(tree, output):
    directory = os.path.dirname(os.path.abspath(output)) or "."
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".solution_", suffix=".xml", dir=directory)
    os.close(fd)
    try:
        tree.write(tmp, encoding="utf-8", xml_declaration=False)
        os.replace(tmp, output)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def verify_output(output, expected_ids, qdec_group):
    tree = ET.parse(output)
    root = tree.getroot()
    element_list = get_element_list(root)

    actual_ids = []
    state_sets = []
    for s in element_list.findall("Set"):
        if not is_state_set(s):
            continue
        if s.get("elementType", "") not in SOLID_ELEMENT_TYPES:
            continue
        state_sets.append(s)
        for e in s.findall("E"):
            actual_ids.append(parse_element_record(e)[0])

    if sorted(actual_ids) != sorted(expected_ids):
        raise RuntimeError(
            "FINAL XML CHECK FAILED: generated state sets do not contain "
            "exactly the physical element IDs."
        )

    loads = root.findall("./Loads/ThermalLoadList/ThermalLoad")
    qloads = [ld for ld in loads if ld.get("uname", "").startswith(QDEC_PREFIX)]
    if len(qloads) != len(expected_ids):
        raise RuntimeError(
            f"FINAL XML CHECK FAILED: expected {len(expected_ids)} Q_dec loads, "
            f"found {len(qloads)}."
        )

    wrong_group = []
    for ld in qloads:
        for sel in ld.findall("Selection"):
            if sel.get("groupName") != qdec_group:
                wrong_group.append(ld.get("uname", ""))
    if wrong_group:
        raise RuntimeError(
            "FINAL XML CHECK FAILED: generated Q_dec loads do not use the "
            f"detected group '{qdec_group}'."
        )

    materials = [
        m for m in root.findall("./MaterialList/Material")
        if is_generated_material(m)
    ]
    if len(materials) != N_STATES:
        raise RuntimeError(
            f"FINAL XML CHECK FAILED: expected {N_STATES} generated ORTHO "
            f"materials, found {len(materials)}."
        )

    # Confirm the detected collector still exists in SetList.
    set_inventory = build_setlist_inventory(root)
    if qdec_group not in set_inventory:
        raise RuntimeError(
            f"FINAL XML CHECK FAILED: Q_dec group '{qdec_group}' no longer exists."
        )

    return len(actual_ids), len(qloads), len(materials)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="TMG/NX 11 T700/M21 ORTHO thermochemical coupling"
    )
    parser.add_argument("--tempf", required=True, help="ASCII TEMPF produced by TMG")
    parser.add_argument("--xml", required=True, help="Input TMG solution XML")
    parser.add_argument("--output", required=True, help="Next TMG solution XML")
    parser.add_argument(
        "--state",
        default="decomposition_state.npz",
        help="Persistent alpha1/alpha2 state file",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=DEFAULT_DT,
        help="Coupling interval [s]",
    )
    args = parser.parse_args()

    if args.dt <= 0.0:
        raise RuntimeError("--dt must be > 0.")

    tree = ET.parse(args.xml)
    root = tree.getroot()

    # --------------------------------------------------------
    # Read all original physical solid elements first.
    # --------------------------------------------------------
    element_list = get_element_list(root)
    all_solid_sets = solid_element_sets(element_list)
    if not all_solid_sets:
        raise RuntimeError("No 3-D solid element sets found in XML.")

    all_element_data, all_element_types, all_element_names = collect_from_element_sets(
        all_solid_sets
    )

    # --------------------------------------------------------
    # Auto-detect the actual thermochemical physical domain.
    # --------------------------------------------------------
    element_data, element_types, element_set_names, source_mode = detect_physical_domain(
        root, all_element_data
    )
    element_ids = sorted(element_data)

    # --------------------------------------------------------
    # Auto-detect the valid NX SetList group used by Q_dec.
    # --------------------------------------------------------
    qdec_group = ensure_qdec_group(root, element_ids)

    print("=" * 76)
    print("TMG / NX 11 T700/M21 ORTHO THERMOCHEMICAL COUPLING")
    print("=" * 76)
    print(f"Input XML       : {args.xml}")
    print(f"TEMPF           : {args.tempf}")
    print(f"Output XML      : {args.output}")
    print(f"State file      : {args.state}")
    print(f"dt              : {args.dt:.9g} s")
    print(f"Physical domain : {source_mode}")
    print(f"Solid elements  : {len(element_ids)}")
    print(f"Element IDs     : {element_ids[0]} .. {element_ids[-1]}")
    print(f"Q_dec group     : {qdec_group}")
    print(f"A1 / Ea1 / n1   : {A1:.7g} 1/s / {Ea1/1000:.7g} kJ/mol / {n1}")
    print(f"A2 / Ea2 / n2   : {A2:.7g} 1/s / {Ea2/1000:.7g} kJ/mol / {n2}")
    print(f"Kcat2           : {KCAT2:.7g}")
    print(f"Qd1 / Qd2       : {Q_DEC1/1000:.5g} / {Q_DEC2/1000:.5g} kJ/kg")
    print(f"epsilon_v/e      : {EPSILON_V:.6f} / {EPSILON_E:.6f}")

    # --------------------------------------------------------
    # TEMPF -> actual physical EIDs.
    # --------------------------------------------------------
    temp_values = parse_temperature_pairs(args.tempf)
    T_C = map_temperatures_to_elements(temp_values, element_data)
    T_K = T_C + 273.15
    print(
        f"TEMPF temp [C]  : min={T_C.min():.4f}, "
        f"max={T_C.max():.4f}, mean={T_C.mean():.4f}"
    )

    # --------------------------------------------------------
    # Advance chemistry.
    # --------------------------------------------------------
    a1, a2 = load_state(args.state, element_ids)
    a1_new, a2_new = rk4_step(T_K, a1, a2, args.dt)

    alpha = decomposition_degree(a1_new, a2_new)

    # Surface radiation coupling: alpha -> effective emissivity -> NX
    # SERadiation property. This is done from the newly advanced chemistry.
    radiation_updates = update_effective_emissivity(root, element_ids, alpha)

    rho = density_from_reactions(a1_new, a2_new)

    # Average rates over this complete coupling interval.  This avoids
    # evaluating Q_dec from only the end-of-step state.
    r1_avg = np.maximum((a1_new - a1) / args.dt, 0.0)
    r2_avg = np.maximum((a2_new - a2) / args.dt, 0.0)
    qdec = q_decomposition_from_average_rates(T_C, r1_avg, r2_avg)

    states_array = np.ceil(alpha / STATE_DA).astype(int)
    states_array = np.clip(states_array, 1, N_STATES)
    states_array[alpha <= 0.0] = 1

    print(
        f"alpha           : min={alpha.min():.6e}, "
        f"max={alpha.max():.6e}, mean={alpha.mean():.6e}"
    )
    print(
        f"density [kg/m3] : min={rho.min():.4f}, max={rho.max():.4f}"
    )
    print(
        f"Q_dec [W/m3]    : min={qdec.min():.6e}, "
        f"max={qdec.max():.6e}, mean={qdec.mean():.6e}"
    )
    print(f"states used     : {sorted(set(states_array.tolist()))}")

    # --------------------------------------------------------
    # Build XML material/state/load data.
    # --------------------------------------------------------
    material_uids = build_material_states(root)
    state_by_eid = {
        eid: int(state)
        for eid, state in zip(element_ids, states_array)
    }

    rebuild_state_sets(
        root,
        element_data,
        element_types,
        state_by_eid,
        material_uids,
    )

    update_qdec_loads(root, element_ids, qdec, qdec_group)

    # Persist only after all XML construction has succeeded.
    save_state_atomic(args.state, element_ids, a1_new, a2_new)

    indent_xml(root)
    write_xml_atomic(tree, args.output)

    # --------------------------------------------------------
    # Final verification.
    # --------------------------------------------------------
    n_state_elems, n_qdec, n_mats = verify_output(
        args.output, element_ids, qdec_group
    )

    print(f"Generated materials : {n_mats}")
    print(f"Generated Q_dec     : {n_qdec}")
    print(f"State-set elements  : {n_state_elems}")
    print(f"Q_dec group kept    : {qdec_group}")
    print(f"Written             : {args.output}")
    print("XML verification    : PASS")
    print("=" * 76)


if __name__ == "__main__":
    main()

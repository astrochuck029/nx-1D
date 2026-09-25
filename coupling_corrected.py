import numpy as np
import argparse
import copy
import os
import re
import tempfile
import xml.etree.ElementTree as ET

R_GAS  = 8.314
DEFAULT_DT = 0.4
N_STATES = 1000
STATE_DA = 0.001

RHO_V = 1575
RHO_E = 1165
DRHO_TOTAL = RHO_V - RHO_E
EPSILON_V = 0.95
EPSILON_E = 0.99

A1 = 7.2193 
A2 = 5.3309e8
Ea1 = 58.7107e3
Ea2 = 146.6740e3
n1 = 1.1125
n2 = 2.0825
KCAT2 = 7.4886
F1 = 0.5540
F2 = 0.9479

Q_DEC1 = 259.54e3
Q_DEC2 = -152.22e3

T_TABLE_C = np.arange(20, 1001, 5)

#Defining the state and material name
SOLID_ELEMENT_TYPES = {"HEXA8", "HEXA20", "TET4", "TET10", "PENTA6", "PENTA15",}
STATE_PREFIX = "3d_mesh_state_"
MAT_PREFIX = "Alpha_State_"
QDEC_PREFIX = "Q_dec_E"
AUTO_GROUP_NAME = "ThermoCouplingDomain"
SERAD_PREFIX = "SERad_E"


#Setting limit for uID (numbers)
MAX_TMG_UID = 99999
MATERIAL_UID_START = 100
QDEC_UID_START = 10000
SERAD_UID_START = 20000

def Cp_virgin(T_C):
    return 2.8773 * T_C + 687.31


def Cp_degraded(T_C):
    return (5.1327e-7 * T_C**3 - 2.0761e-3 * T_C**2 + 2.599 * T_C + 662.53)


def k_virgin_inplane(T_C):
    return 7.4675e-3 * T_C + 2.7811


def k_virgin_thru(T_C):
    return 1.1113e-3 * T_C + 0.61391


def k_degraded_inplane(T_C):
    return -3.5481e-6 * T_C**2 + 6.4898e-3 * T_C + 2.3005


def k_degraded_thru(T_C):
    return (7.4228e-10 * T_C**3 - 4.1903e-7 * T_C**2 + 2.3397e-4 * T_C + 7.7211e-2)


def mixture_properties(T_C, alpha):
    alpha = float(np.clip(alpha, 0.0, 1.0))
    T_C = np.asarray(T_C, dtype=float)

    rho = RHO_V - DRHO_TOTAL * alpha
    cp = (1 - alpha) * Cp_virgin(T_C) + alpha * Cp_degraded(T_C)
    kip = (1 - alpha) * k_virgin_inplane(T_C) + alpha * k_degraded_inplane(T_C) #Kinplane update
    kth = (1 - alpha) * k_virgin_thru(T_C) + alpha * k_degraded_thru(T_C) #kthru_plane update
    return rho, cp, kip, kip.copy(), kth

#Reaction Kinetics and Decompostion
#Gas enthalpy rise 
def gas_enthalpy_rise(T_C, T0_C=30): #T_init = T0_C hardcoded to match with the NX initialization temperature
    T_C = np.asarray(T_C, dtype=float)

    def H(T):
        return ((3.5977e-7 / 4) * T**4 - (9.2485e-4 / 3) * T**3 + (1.0610 / 2) * T**2 + 1256.6 * T )

    return np.maximum(H(T_C) - H(T0_C), 0)
#Two step Arrehenous Reaction
def reaction_rates(T_K, a1, a2):
    T_K = np.maximum(np.asarray(T_K, dtype=float), 1)
    a1 = np.asarray(a1, dtype=float)
    a2 = np.asarray(a2, dtype=float)
    remaining = np.clip(1-a1-a2, 0, 1)
    r1 = A1 * remaining**n1 * np.exp(-Ea1 / (R_GAS * T_K))
    r2 = (A2 * remaining**n2 * np.exp(-Ea2 / (R_GAS * T_K)) * (1 + KCAT2 * np.clip(a2, 0, 1)))
    return r1, r2

def rk4_step(T_K, a1, a2, dt):
    if dt<=0:
        raise ValueError("!!!Keep dt >0")

    k1a, k1b = reaction_rates(T_K, a1, a2)
    k2a, k2b = reaction_rates(T_K, a1+(0.5 * k1a * dt), a2+(0.5 * k1b * dt))
    k3a, k3b = reaction_rates(T_K, a1+(0.5 * k2a * dt), a2+(0.5 * k2b * dt))
    k4a, k4b = reaction_rates(T_K, a1+(dt * k3a), a2 + (dt* k3b))

    new_a1 = a1 + (dt/6)*(k1a + 2*k2a + 2*k3a + k4a)
    new_a2 = a2 + (dt/6)*(k1b + 2*k2b + 2*k3b + k4b)

    new_a1 = np.clip(new_a1, 0, 1)
    new_a2 = np.clip(new_a2,0, 1)

    total = new_a1 + new_a2
    over  = total>1
    if np.any(over):
        scale = np.ones_like(total)
        scale[over] = 1/total[over]
        new_a1 *= scale
        new_a2 *= scale
    return new_a1, new_a2

def decomposition_degree(a1, a2):
    return np.clip(F1*a1 + F2*a2, 0, 1)

def density_from_reaction(a1, a2):
    rho = RHO_V - DRHO_TOTAL*decomposition_degree(a1, a2)
    return np.clip(rho, RHO_E, RHO_V)

def solid_enthalpy_rise(T_C, a1, a2, T0_C=30):
    def hs_v(T): 
        return 687.31 * T + (2.8773 / 2.0) * T**2
    
    def hs_c(T): 
        return (662.53 * T + (2.599 / 2) * T**2 - (2.0761e-3 / 3) * T**3 + (5.132e-7 / 4) * T**4)
    #alpha = F1 * a1 + F2 * a2
    hs = (1 - decomposition_degree(a1, a2)) * (hs_v(T_C) - hs_v(T0_C)) + decomposition_degree(a1, a2) * (hs_c(T_C) - hs_c(T0_C))
    return np.maximum(hs, 0)

def q_dec_from_avg_rates(T_C, r1_avg, r2_avg, a1, a2):
    q_chem = RHO_V * (F1 * Q_DEC1 * r1_avg + F2 * Q_DEC2 * r2_avg)
    
    drho1 = F1 * DRHO_TOTAL
    drho2 = F2 * DRHO_TOTAL
    drho_g = drho1 * r1_avg + drho2 * r2_avg
    
    dh_gas = gas_enthalpy_rise(T_C, T0_C=30)
    dh_solid = solid_enthalpy_rise(T_C, a1, a2, T0_C=30)

    return q_chem - drho_g * (dh_solid - dh_gas)

# def q_dec_from_avg_rates(T_C, r1_avg, r2_avg):
#     drho1 = F1 * DRHO_TOTAL
#     drho2 = F2 * DRHO_TOTAL

#     q_pyro =  Q_DEC1 * drho1 * r1_avg + Q_DEC2 * drho2 * r2_avg
#     drho_g = drho1 * r1_avg + drho2 * r2_avg
#     dh_gas = gas_enthalpy_rise(T_C, T0_C=30)

#     return q_pyro - drho_g* dh_gas


"""Some Xml Helper Function"""
# Now I have to crate some xml file helper functions
#This fuctions creats an specific format [scientic notation] to match with nx xml
def fmt(value):
    return f" {float(value):.7E}"

def get_element_list(root):
    node = root.find("ElementList")
    if node is None:
        raise RuntimeError("Could't find <ElementList>")
    return node

def get_set_list(root):
    node = root.find("SetList")
    if node is None:
        raise RuntimeError("Could not find <SetList>")
    return node

def parse_element_record(e):
    text = (e.text or "").strip()
    if not text:
        raise RuntimeError("There is <E> entry with no element data")

    values = text.split()

    try:
        eid = int(values[0])
    except ValueError as exc:
        raise RuntimeError(f"Invalide EID in <E>:{text!r}") from exc

    if eid <=0:
        raise RuntimeError(f"Not a valid EID {eid}")

    return eid, values

def is_state_set(s):
    return s.get("uname", "").startswith(STATE_PREFIX)

def is_generated_material(m):
    return m.get("uname", "").startswith(MAT_PREFIX)

def is_generated_qdec(load):
    uname = load.get("uname", "")
    return uname == "Q_dec" or uname.startswith(QDEC_PREFIX)


def solid_element_sets(element_list):
    return [
        s for s in element_list.findall("Set")
        if s.get("elementType", "") in SOLID_ELEMENT_TYPES
        and s.findall("E")
    ]

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
                    f"Duplicate solid element ID {eid} found in XML "
                    f"(set '{names[eid]}' and '{uname}')."
                )

            data[eid] = values
            types[eid] = etype
            names[eid] = uname

    if not data:
        raise RuntimeError("No solid elements found in XML.")

    return data, types, names

def setlist_selection_ids(s):
    ids = []

    for sel in s.findall("Selection"):
        for item in sel.findall("_"):
            text = (item.text or "").strip()
            if text:
                try:
                    ids.append(int(text))
                except ValueError:
                    pass

        for item in sel.findall("el"):
            text = (item.text or "").strip()
            if text:
                try:
                    ids.append(int(text))
                except ValueError:
                    pass

    return sorted(set(ids))


def build_setlist_collector(root):
    collector =  {}
    for s in get_set_list(root).findall("Set"):
        name = s.get("uname", "")
        if name:
            collector[name] = set(setlist_selection_ids(s))

    return collector

#The above fuction build_setlist_collector resturn this output for current xml
# 'Heated_Faces': {1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241},

def thermal_load_group_names(root):
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

def group_name_from_existing_qdec(root):
    loads = root.find("./Loads/ThermalLoadList")

    if loads is None:
        return None

    for load in loads.findall("ThermalLoad"):
        if load.get("uname", "") != "Q_dec":   #Q_dec Hardcoded here
            continue

        for sel in load.findall("Selection"):
            name = sel.get("groupName")
            if name:
                return name
    return None

def find_physical_domain(root):
    element_list = get_element_list(root)
    all_solid_sets = solid_element_sets(element_list)

    if not all_solid_sets:
        raise RuntimeError("No 3-D solid element sets found in XML.")

    state_sets = [
        s for s in all_solid_sets
        if is_state_set(s)
    ]

    if state_sets:
        data, types, names = collect_from_element_sets(state_sets)
        print("Physical domain source : existing generated state sets")
        return data, types, names, "generated state sets"

    # Otherwise use a complete SetList collector if one exists.
    solid_ids = set()

    for s in all_solid_sets:
        for e in s.findall("E"):
            solid_ids.add(parse_element_record(e)[0])

    set_inventory = build_setlist_collector(root)
    load_groups = thermal_load_group_names(root)
    q_group = group_name_from_existing_qdec(root)

    candidates = []

    for name, ids in set_inventory.items():
        overlap = ids & solid_ids
        if not overlap:
            continue

        candidates.append(
            (
                len(overlap) == len(solid_ids),
                q_group == name,
                name in load_groups,
                len(overlap),
                name,
                overlap,
            )
        )

    if candidates:
        candidates.sort(
            reverse=True,
            key=lambda x: (x[0], x[1], x[2], x[3], x[4])
        )

        _, is_q, referenced, _, name, overlap = candidates[0]

        selected_sets = [
            s for s in all_solid_sets
            if any(
                parse_element_record(e)[0] in overlap
                for e in s.findall("E")
            )
        ]

        data, types, names = collect_from_element_sets_filtered(
            selected_sets, overlap
        )

        reason = "SetList collector"
        if is_q:
            reason += " referenced by existing Q_dec"
        elif referenced:
            reason += " referenced by thermal load"

        print(f"Physical group        : {name}")
        print(f"Physical group source : {reason}")

        return data, types, names, f"collector '{name}'"

    # No complete collector: use all solid elements.
    data, types, names = collect_from_element_sets(all_solid_sets)
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
                    f"Duplicate physical solid element ID {eid} "
                    f"while using a collector."
                )

            data[eid] = values
            types[eid] = etype
            names[eid] = uname

    if set(data) != allowed_ids:
        missing = sorted(allowed_ids - set(data))
        preview = ", ".join(map(str, missing[:20]))

        if len(missing) > 20:
            preview += ", ..."

        raise RuntimeError(f"Collector references solid element IDs that cannot be found in ElementList: {preview}")

    return data, types, names

def ensure_qdec_group(root, element_ids):
    element_ids = set(element_ids)
    inventory = build_setlist_collector(root)

    existing_q = group_name_from_existing_qdec(root)

    if existing_q and existing_q in inventory:
        if element_ids.issubset(inventory[existing_q]):
            print(f"Q_dec group            : {existing_q}")
            return existing_q

    full = [
        name for name, ids in inventory.items()
        if element_ids.issubset(ids)
    ]

    if full:
        refs = set(thermal_load_group_names(root))
        full.sort(key=lambda n: (n not in refs, len(inventory[n]), n))

        name = full[0]
        print(f"Q_dec group            : {name}")
        return name

    ranked = []

    refs = set(thermal_load_group_names(root))

    for name, ids in inventory.items():
        overlap = len(ids & element_ids)

        if overlap:
            ranked.append(
                (overlap == len(element_ids), overlap, name in refs, name)
            )

    if ranked:
        ranked.sort(reverse=True)

        name = ranked[0][3]
        print(f"Q_dec group            : {name}")
        return name

    set_list = get_set_list(root)

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

    s = ET.SubElement(
        set_list,
        "Set",
        {"uid": str(uid), "uname": name}
    )

    sel = ET.SubElement(
        s,
        "Selection",
        {"step": "1", "class": "Elements"}
    )

    for eid in sorted(element_ids):
        ET.SubElement(sel, "_").text = str(eid)

    print(f"Q_dec group            : {name} (created)")
    return name
######################################################################################################
"""Material Unique identifiers (UID) using xml.etree API helper Functions"""
#I need to look for all uid and map with the element naterial properties
def all_uids(root, exclude_generated=False):
    used = set()
    for elem in root.iter():
        if exclude_generated:
            if elem.tag == "Material" and is_generated_material(elem):
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
    return max(all_uids(root), default=0)

def allocate_uid_block(root, count, preferred_start):
    if count <= 0:
        raise ValueError("count must be > 0")

    used = all_uids(root, exclude_generated=True)
    start = max(1, int(preferred_start))

    for candidate in (start, max_uid(root) + 1):
        candidate = max(1, candidate)
        end = candidate + count - 1

        if end > MAX_TMG_UID:
            continue

        block = range(candidate, end + 1)

        if all(uid not in used for uid in block):
            return candidate

    raise RuntimeError(
        f"Could not define {count} unique TMG UIDs below {MAX_TMG_UID}."
    )

#########################################################################
#Material Properties defination in xml

def replace_property_with_table(material, name, values):
    old = None

    for p in material.findall("Property"):
        if p.get("name") == name:
            old = p
            break

    new = ET.Element(
        "Property",
        {"name": name, "type": "XYTable", "x": "Temperature"}
    )

    for T, value in zip(T_TABLE_C, values):
        ET.SubElement(new, "_").text = (
            f"{T: .7E}    {float(value): .7E}"
        )

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
        old = ET.SubElement( material, "Property",
            {"name": name, "type": "Constant"}
        )
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

    candidates = [
        m for m in material_list.findall("Material")
        if m.get("type", "").upper() == "ORTHO"
        and not is_generated_material(m)
    ]

    if not candidates:
        raise RuntimeError("No original ORTHO material found in current xml")

    return candidates[0]


def build_material_states(root):
    material_list = root.find("MaterialList")

    if material_list is None:
        raise RuntimeError("Could not find <MaterialList> in XML.")

    template = find_ortho_template(root)

    for m in list(material_list.findall("Material")):
        if is_generated_material(m):
            material_list.remove(m)

    base_uid = allocate_uid_block(root, N_STATES, MATERIAL_UID_START)

    state_uids = {}

    for state in range(1, N_STATES + 1):
        alpha = state * STATE_DA

        rho, cp, kx, ky, kz = mixture_properties(T_TABLE_C, alpha)

        mat = copy.deepcopy(template)

        mat.set("uid", str(base_uid + state - 1))
        mat.set("uname", f"{MAT_PREFIX}{state:02d}")
        mat.set("type", "ORTHO")

        replace_constant_property(mat, "Mass Density", rho)

        replace_property_with_table( mat, "Specific Heat", cp)

        replace_property_with_table( mat, "Thermal Conductivity X", kx )

        replace_property_with_table( mat, "Thermal Conductivity Y", ky)

        replace_property_with_table(mat, "Thermal Conductivity Z", kz )

        material_list.append(mat)
        state_uids[state] = base_uid + state - 1

    return state_uids

#####################################################################
# State set Replacement Fucntions       

def rebuild_state_sets(
    root,
    element_data,
    element_types,
    states,
    material_uids,
):
    element_list = get_element_list(root)
    expected = sorted(element_data)

    physical_ids = set(expected)
    for s in list(element_list.findall("Set")):
        if s.get("elementType", "") not in SOLID_ELEMENT_TYPES:
            continue

        ids_in_set = {
            parse_element_record(e)[0]
            for e in s.findall("E")
        }

        if ids_in_set & physical_ids:
            element_list.remove(s)

    for state in range(1, N_STATES + 1):
        by_type = {}

        for eid in expected:
            if int(states[eid]) != state:
                continue

            by_type.setdefault(
                element_types[eid],
                []
            ).append(eid)

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

    seen = []

    for s in element_list.findall("Set"):
        if not is_state_set(s):
            continue

        if s.get("elementType", "") not in SOLID_ELEMENT_TYPES:
            continue

        for e in s.findall("E"):
            seen.append(parse_element_record(e)[0])

    if len(seen) != len(set(seen)):
        raise RuntimeError(
            "State-set rebuild failed: duplicate physical element IDs."
        )

    if sorted(seen) != expected:
        raise RuntimeError(
            "State-set rebuild failed: physical solid EIDs changed."
        )

def create_qdec_load(uid, element_id, q_value, group_name):
    load =  ET.Element("ThermalLoad",
                       {
                           "uid": str(uid),
                           "uname": f"{QDEC_PREFIX}{element_id}",
                           "type" : "Heat Generation"
                       },)
    ET.SubElement(load, "Description")

    p = ET.SubElement(
        load,
        "Property",
        {"name": "Heat Generation"}
    )

    ET.SubElement(p, "Value").text = fmt(q_value)

    for name, value in [
        ("Control Heater", "0"),
        ("Thermostat", "-1"),
        ("Specify Layer to Apply to", "0"),
        ("Apply to", "0"),
        ("Layer Number", "1"),
    ]:
        p = ET.SubElement(
            load,
            "Property",
            {"name": name}
        )
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
        raise RuntimeError("Could not find <Loads>/<ThermalLoadList>")
    
    for load in list(loads.findall("ThermalLoad")):
        if is_generated_qdec(load):
            loads.remove(load)

    uid = allocate_uid_block(
        root,
        len(element_ids),
        QDEC_UID_START
    )

    for eid, q in zip(element_ids, qdec):
        loads.append(create_qdec_load(uid, eid, q, group_name))
        uid += 1


def is_generated_serad(rad):
    uname = rad.get("uname", "")
    return (
        uname.startswith("SERad_F_E")
        or uname.startswith("SERad_R_E")
        or uname.startswith(SERAD_PREFIX)
    )


def _selection_data(rad):
    group_name = None
    eids = []

    for sel in rad.findall("Selection"):
        group_name = sel.get("groupName") or group_name

        for item in sel.findall("el"):
            text = (item.text or "").strip()
            if text:
                eids.append(int(text))

        for item in sel.findall("_"):
            text = (item.text or "").strip()
            if text:
                eids.append(int(text))

    return group_name, sorted(set(eids))


def _find_serad_by_names(root, names):
    rad_list = root.find("./Constraints/SERadiationList")
    if rad_list is None:
        return None, [], None

    names = set(names)

    for rad in rad_list.findall("SERadiation"):
        if rad.get("uname", "") in names:
            group_name, eids = _selection_data(rad)
            if group_name and eids:
                return group_name, eids, rad

    return None, [], None


def _find_generated_serad(root, prefix):
    rad_list = root.find("./Constraints/SERadiationList")
    if rad_list is None:
        return None, [], None

    group_name = None
    eids = set()
    template = None

    for rad in rad_list.findall("SERadiation"):
        uname = rad.get("uname", "")
        if not uname.startswith(prefix):
            continue

        if template is None:
            template = rad

        group, ids = _selection_data(rad)
        if group:
            group_name = group
        eids.update(ids)

    if group_name and eids and template is not None:
        return group_name, sorted(eids), template

    return None, [], None


def find_serad_domain(root, default_uname="Radiation_front"):
    """Find one radiation surface and its template load.

    The original TMG radiation load is used as the template so that its
    environment temperature, view factor, temperature type, etc. are kept
    unchanged.  After the first coupling step the original load is replaced
    by per-element generated loads, so those are also recognized here.
    """
    source_names = [default_uname]

    if default_uname == "Radiation_front":
        source_names += ["Radiation_heated", "Radiation_heated_face"]
        generated_prefix = "SERad_F_E"
    else:
        source_names += ["Radiation_back", "Radiation_back_face", "Radiation_rear_face"]
        generated_prefix = "SERad_R_E"

    group_name, eids, template = _find_serad_by_names(root, source_names)
    if group_name and eids:
        return group_name, eids, template

    return _find_generated_serad(root, generated_prefix)


def _set_property_value(load, property_name, value):
    for prop in load.findall("Property"):
        if prop.get("name") == property_name:
            value_node = prop.find("Value")
            if value_node is None:
                value_node = ET.SubElement(prop, "Value")
            value_node.text = fmt(value)
            return

    prop = ET.SubElement(load, "Property", {"name": property_name})
    ET.SubElement(prop, "Value").text = fmt(value)


def create_serad_load(uid, element_id, ep_value, group_name, template, prefix):
    """Create one per-element radiation injection from the original TMG load."""
    if template is None:
        raise RuntimeError(
            f"No SERadiation template was found for group '{group_name}'."
        )

    rad = copy.deepcopy(template)
    rad.set("uid", str(uid))
    rad.set("uname", f"{prefix}{element_id}")
    rad.set("type", "Simple Radiation to Environment")

    # Keep every original radiation property and change only the effective
    # emissivity that is being coupled from the decomposition model.
    _set_property_value(rad, "Effective Emissivity", ep_value)

    # Replace the original selection by exactly one element, just like Q_dec.
    for sel in list(rad.findall("Selection")):
        rad.remove(sel)

    sel = ET.SubElement(
        rad,
        "Selection",
        {"step": "1", "groupName": group_name}
    )
    ET.SubElement(sel, "el").text = str(element_id)

    return rad


def update_serad_loads(
    root,
    element_ids,
    ep_values,
    group_name,
    template,
    prefix,
    source_uname=None,
):
    """Inject one dynamic SERadiation load per surface element.

    This follows the same pattern as Q_dec:
      element 1 -> SERad_F_E1
      element 2 -> SERad_F_E2
      ...

    The same is done independently for the rear surface with SERad_R_E*.
    """
    constraints = root.find("Constraints")
    if constraints is None:
        constraints = ET.SubElement(root, "Constraints")

    rad_list = constraints.find("SERadiationList")
    if rad_list is None:
        rad_list = ET.SubElement(constraints, "SERadiationList")

    # Remove only previously generated coupling loads.
    for rad in list(rad_list.findall("SERadiation")):
        uname = rad.get("uname", "")
        if uname.startswith(prefix) or (
            prefix == "SERad_F_E" and uname.startswith(SERAD_PREFIX)
        ):
            rad_list.remove(rad)

    # Remove the original source load only after its properties have been
    # captured in 'template'.  This prevents the fire/ambient temperature or
    # view-factor settings from being replaced by hard-coded values.
    if source_uname:
        for rad in list(rad_list.findall("SERadiation")):
            if rad.get("uname", "") == source_uname:
                rad_list.remove(rad)

    uid = allocate_uid_block(
        root,
        len(element_ids),
        SERAD_UID_START
    )

    for eid, ep in zip(element_ids, ep_values):
        rad_list.append(
            create_serad_load(
                uid,
                eid,
                ep,
                group_name,
                template,
                prefix,
            )
        )
        uid += 1

#TEMPF Operation
_NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")

def parse_temperature_pairs(filename):
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
        raise RuntimeError("No tempr record found in tempf file")

    return values

def map_temperatures_to_elements(temp_values, element_data):
    eids = sorted(element_data)
    missing = [
        eid for eid in eids
        if eid not in temp_values
    ]

    if missing:
        preview = ", ".join(map(str, missing[:30]))

        if len(missing) > 30:
            preview += ", ..."

        raise RuntimeError(
            "TEMPF/XML ELEMENT-ID MISMATCH\n"
            f"  XML physical elements : {len(eids)}\n"
            f"  TEMPF matched         : {len(eids) - len(missing)}\n"
            f"  Missing IDs           : {preview}"
        )
    return np.array([temp_values[eid] for eid in eids], dtype = float)

#Updating the states
def load_state(filename, element_ids):
    n = len(element_ids)

    if not os.path.exists(filename):
        print(f"No previous alpha state. Initializing with alpha=0 for {n} elements.")

        return np.zeros(n), np.zeros(n)

    with np.load(filename, allow_pickle=False) as data:
        required = {"element_ids", "a1", "a2"}

        if not required.issubset(data.files):
            raise RuntimeError(f"Incompatible state file: {filename}")

        old_ids = np.asarray(data["element_ids"], dtype=np.int64)

        a1 = np.asarray(data["a1"], dtype=float)

        a2 = np.asarray(data["a2"], dtype=float)

    ids = np.asarray(element_ids, dtype=np.int64)

    if old_ids.shape != ids.shape or not np.array_equal(old_ids, ids):
        raise RuntimeError("Saved decomposition state does not match the current eID")

    if a1.shape != ids.shape or a2.shape != ids.shape:
        raise RuntimeError("Saved alpha arrays do not match the current mesh")

    a1 = np.clip(a1, 0.0, 1.0)
    a2 = np.clip(a2, 0.0, 1.0)

    total = a1 + a2
    over = total > 1.0

    if np.any(over):
        a1[over] /= total[over]
        a2[over] /= total[over]

    return a1, a2

def save_state_atomic(filename, element_ids, a1, a2):
    directory = os.path.dirname(
        os.path.abspath(filename)
    ) or "."

    os.makedirs(directory, exist_ok=True)

    fd, tmp = tempfile.mkstemp(
        prefix=".decomp_",
        suffix=".npz",
        dir=directory
    )

    os.close(fd)

    try:
        with open(tmp, "wb") as f:  # Use a file handler so numpy does not append another .npz
            np.savez(
                f,
                element_ids=np.asarray(
                    element_ids,
                    dtype=np.int64
                ),
                a1=np.asarray(a1, dtype=float),
                a2=np.asarray(a2, dtype=float),
            )

        os.replace(tmp, filename)

    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

##########################################################################################################
#XML verification
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

    fd, tmp = tempfile.mkstemp(
        prefix=".solution_",
        suffix=".xml",
        dir=directory
    )

    os.close(fd)

    try:
        tree.write(tmp, encoding="utf-8", xml_declaration=False)

        os.replace(tmp, output)

    finally:
        if os.path.exists(tmp):
            os.remove(tmp)



def verify_output(output, expected_ids, qdec_group, expected_front_ids=None, expected_rear_ids=None, front_group=None, rear_group=None):
    tree = ET.parse(output)
    root = tree.getroot()

    element_list = get_element_list(root)
    expected = sorted(expected_ids)
    physical = set(expected)

    state_ids = []

    for s in element_list.findall("Set"):
        if not is_state_set(s):
            continue

        if s.get("elementType", "") not in SOLID_ELEMENT_TYPES:
            continue

        for e in s.findall("E"):
            state_ids.append(
                parse_element_record(e)[0]
            )

    if len(state_ids) != len(set(state_ids)):
        raise RuntimeError("FINAL XML CHECK FAILED: duplicate EIDs inside generated state sets")

    if sorted(state_ids) != expected:
        raise RuntimeError("FINAL XML CHECK FAILED: generated state sets do not contain exactly the physical element IDs")

    duplicate_old = []

    for s in element_list.findall("Set"):
        if is_state_set(s):
            continue

        if s.get("elementType", "") not in SOLID_ELEMENT_TYPES:
            continue

        for e in s.findall("E"):
            eid = parse_element_record(e)[0]

            if eid in physical:
                duplicate_old.append(
                    (eid, s.get("uname", ""))
                )

    if duplicate_old:
        preview = ", ".join(
            f"{eid}:{name}"
            for eid, name in duplicate_old[:20]
        )

        raise RuntimeError(f"FINAL XML CHECK FAILED: physical EIDs remain in old solid ElementList sets: {preview}")

    loads = root.findall(
        "./Loads/ThermalLoadList/ThermalLoad"
    )

    qloads = [
        ld for ld in loads
        if ld.get("uname", "").startswith(QDEC_PREFIX)
    ]

    if len(qloads) != len(expected):
        raise RuntimeError(
            f"FINAL XML CHECK FAILED: expected {len(expected)} "
            f"generated Q_dec loads, found {len(qloads)}."
        )

    qids = []

    for ld in qloads:
        uname = ld.get("uname", "")
        match = re.fullmatch(
            rf"{re.escape(QDEC_PREFIX)}(\d+)",
            uname
        )

        if not match:
            raise RuntimeError(f"FINAL XML CHECK FAILED: invalid Q_dec name {uname}")

        qids.append(int(match.group(1)))

        for sel in ld.findall("Selection"):
            if sel.get("groupName") != qdec_group:
                raise RuntimeError(
                    f"FINAL XML CHECK FAILED: {uname} uses the wrong "
                    "Q_dec group."
                )

    if sorted(qids) != expected:
        raise RuntimeError(
            "FINAL XML CHECK FAILED: generated Q_dec loads do not "
            "cover exactly the physical EIDs."
        )

    if any(ld.get("uname", "") == "Q_dec" for ld in loads):
        raise RuntimeError(
            "FINAL XML CHECK FAILED: original Q_dec load remains"
        )

    #Total numebr of generated materials
    materials = [m for m in root.findall("./MaterialList/Material")
        if is_generated_material(m)]

    if len(materials) != N_STATES:
        raise RuntimeError(
            f"FINAL XML CHECK FAILED: expected {N_STATES} "
            f"generated materials, found {len(materials)}."
        )

    if qdec_group not in build_setlist_collector(root):
        raise RuntimeError(
            f"FINAL XML CHECK FAILED: Q_dec group "
            f"'{qdec_group}' no longer exists."
        )

    rad_loads = root.findall("./Constraints/SERadiationList/SERadiation")

    def _generated_for_prefix(prefix):
        return [
            r for r in rad_loads
            if r.get("uname", "").startswith(prefix)
        ]

    front_generated = _generated_for_prefix("SERad_F_E")
    rear_generated = _generated_for_prefix("SERad_R_E")

    if expected_front_ids is not None:
        if len(front_generated) != len(expected_front_ids):
            raise RuntimeError(
                f"FINAL XML CHECK FAILED: expected {len(expected_front_ids)} "
                f"front SERadiation loads, found {len(front_generated)}."
            )

        front_ids = []
        for rad in front_generated:
            group, ids = _selection_data(rad)
            if front_group is not None and group != front_group:
                raise RuntimeError(
                    f"FINAL XML CHECK FAILED: front SERadiation {rad.get('uname')} "
                    "uses the wrong group."
                )
            front_ids.extend(ids)

        if sorted(front_ids) != sorted(expected_front_ids):
            raise RuntimeError(
                "FINAL XML CHECK FAILED: front SERadiation loads do not "
                "cover exactly the front surface EIDs."
            )

    if expected_rear_ids is not None:
        if len(rear_generated) != len(expected_rear_ids):
            raise RuntimeError(
                f"FINAL XML CHECK FAILED: expected {len(expected_rear_ids)} "
                f"rear SERadiation loads, found {len(rear_generated)}."
            )

        rear_ids = []
        for rad in rear_generated:
            group, ids = _selection_data(rad)
            if rear_group is not None and group != rear_group:
                raise RuntimeError(
                    f"FINAL XML CHECK FAILED: rear SERadiation {rad.get('uname')} "
                    "uses the wrong group."
                )
            rear_ids.extend(ids)

        if sorted(rear_ids) != sorted(expected_rear_ids):
            raise RuntimeError(
                "FINAL XML CHECK FAILED: rear SERadiation loads do not "
                "cover exactly the rear surface EIDs."
            )

    if any(
        r.get("uname", "") in {
            "Radiation_front",
            "Radiation_rear",
            "Radiation_back",
            "Radiation_back_face",
            "Radiation_heated",
            "Radiation_heated_face",
            "Radiation_rear_face",
        }
        for r in rad_loads
    ):
        raise RuntimeError(
            "FINAL XML CHECK FAILED: an original front/rear Radiation load remains."
        )

    return (
        len(state_ids),
        len(qloads),
        len(materials),
        len(front_generated),
        len(rear_generated),
    )


def main():
    parser = argparse.ArgumentParser(
        description="TMG/NX11 Thermochemical Coupling"
    )

    parser.add_argument(
        "--tempf",
        required=True,
        help="ASCII TEMPF produced by TMG"
    )

    parser.add_argument(
        "--xml",
        required=True,
        help="Input TMG solution XML"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Next TMG solution XML"
    )

    parser.add_argument(
        "--state",
        default="decomposition_state.npz",
        help="Persistent decomposition state"
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=DEFAULT_DT,
        help="Coupling interval [s]"
    )

    args = parser.parse_args()

    if args.dt <= 0.0:
        raise RuntimeError("--dt must be > 0.")

    tree = ET.parse(args.xml)
    root = tree.getroot()

    element_data, element_types, _, source_mode = (find_physical_domain(root))

    element_ids = sorted(element_data)

    qdec_group = ensure_qdec_group(root, element_ids)

    print("/" * 76)
    print("TMG / NX 11 T700/M21 THERMOCHEMICAL COUPLING")
    print("/" * 76)
    print(f"Input XML       : {args.xml}")
    print(f"TEMPF           : {args.tempf}")
    print(f"Output XML      : {args.output}")
    print(f"State file      : {args.state}")
    print(f"dt              : {args.dt:.9g} s")
    print(f"Physical domain : {source_mode}")
    print(f"Solid elements  : {len(element_ids)}")
    print(f"Element IDs     : {element_ids[0]} .. {element_ids[-1]}")
    print(f"Q_dec group     : {qdec_group}")

    temp_values = parse_temperature_pairs(args.tempf)

    T_C = map_temperatures_to_elements(temp_values, element_data)

    T_K = T_C + 273.15

    print(
        f"TEMPF temp [C]  : min={T_C.min():.4f}, "
        f"max={T_C.max():.4f}, "
        f"mean={T_C.mean():.4f}"
    )

    a1, a2 = load_state(
        args.state,
        element_ids
    )

    a1_new, a2_new = rk4_step(
        T_K,
        a1,
        a2,
        args.dt
    )

    alpha = decomposition_degree(
        a1_new,
        a2_new
    )

    rho = density_from_reaction(
        a1_new,
        a2_new
    )

    r1_avg = np.maximum(
        (a1_new - a1) / args.dt,
        0.0
    )

    r2_avg = np.maximum(
        (a2_new - a2) / args.dt,
        0.0
    )

    qdec = q_dec_from_avg_rates(T_C, r1_avg, r2_avg, a1, a2)

    states_array = np.ceil(
        alpha / STATE_DA
    ).astype(int)

    states_array = np.clip(
        states_array,
        1,
        N_STATES
    )

    states_array[alpha <= 0.0] = 1

    print(
        f"alpha           : min={alpha.min():.6e}, "
        f"max={alpha.max():.6e}, "
        f"mean={alpha.mean():.6e}"
    )

    print(
        f"density [kg/m3] : min={rho.min():.4f}, "
        f"max={rho.max():.4f}"
    )

    print(
        f"Q_dec [W/m3]    : min={qdec.min():.6e}, "
        f"max={qdec.max():.6e}, "
        f"mean={qdec.mean():.6e}"
    )

    print(
        f"states used     : "
        f"{sorted(set(states_array.tolist()))}"
    )

    # Surface emissivity updates
    # Each surface is injected independently, using one SERadiation load per
    # element, exactly like the existing Q_dec per-element injection.
    front_group, front_eids, front_template = find_serad_domain(
        root,
        default_uname="Radiation_front"
    )

    rear_group, rear_eids, rear_template = find_serad_domain(
        root,
        default_uname="Radiation_rear"
    )

    if not front_eids:
        raise RuntimeError(
            "Could not find the heated/front radiation surface. "
            "Expected Radiation_front or its generated SERad_F_E* loads."
        )

    if not rear_eids:
        raise RuntimeError(
            "Could not find the rear radiation surface. "
            "Expected Radiation_rear (or Radiation_back) or its generated "
            "SERad_R_E* loads."
        )

    print(f"Front SERad group       : {front_group}")
    print(f"Front surface elements : {len(front_eids)} ({front_eids[0]} .. {front_eids[-1]})")
    print(f"Rear SERad group        : {rear_group}")
    print(f"Rear surface elements  : {len(rear_eids)} ({rear_eids[0]} .. {rear_eids[-1]})")

    missing_front = [eid for eid in front_eids if eid not in temp_values]
    if missing_front:
        preview = ", ".join(map(str, missing_front[:20]))
        raise RuntimeError(
            f"TEMPF file missing temperatures for front surface elements: {preview}"
        )

    missing_rear = [eid for eid in rear_eids if eid not in temp_values]
    if missing_rear:
        preview = ", ".join(map(str, missing_rear[:20]))
        raise RuntimeError(
            f"TEMPF file missing temperatures for rear surface elements: {preview}"
        )

    # Keep independent reaction histories for the two surfaces because their
    # temperature histories are generally different.
    base_state, ext_state = os.path.splitext(args.state)
    state_ext = ext_state if ext_state else ".npz"
    front_state_file = f"{base_state}_front{state_ext}"
    rear_state_file = f"{base_state}_rear{state_ext}"

    # ----------------------- Heated / front surface -----------------------
    T_C_front = np.array(
        [temp_values[eid] for eid in front_eids],
        dtype=float
    )
    T_K_front = T_C_front + 273.15

    a1_front, a2_front = load_state(
        front_state_file,
        front_eids
    )

    a1_front_new, a2_front_new = rk4_step(
        T_K_front,
        a1_front,
        a2_front,
        args.dt
    )

    alpha_front = decomposition_degree(
        a1_front_new,
        a2_front_new
    )

    ep_front = np.clip(
        (1.0 - alpha_front) * EPSILON_V
        + alpha_front * EPSILON_E,
        0.0,
        1.0,
    )

    save_state_atomic(
        front_state_file,
        front_eids,
        a1_front_new,
        a2_front_new
    )

    # The original Radiation_front load is used as the template. Therefore
    # its fire/ambient temperature and view-factor settings are preserved.
    update_serad_loads(
        root,
        front_eids,
        ep_front,
        front_group,
        front_template,
        "SERad_F_E",
        source_uname=(
            front_template.get("uname", "")
            if front_template is not None
            and front_template.get("uname", "").startswith("Radiation_")
            else None
        ),
    )

    print(
        f"Front emissivity        : min={ep_front.min():.6f}, "
        f"max={ep_front.max():.6f}, mean={ep_front.mean():.6f}"
    )

    # ------------------------- Rear surface -------------------------------
    T_C_rear = np.array(
        [temp_values[eid] for eid in rear_eids],
        dtype=float
    )
    T_K_rear = T_C_rear + 273.15

    a1_rear, a2_rear = load_state(
        rear_state_file,
        rear_eids
    )

    a1_rear_new, a2_rear_new = rk4_step(
        T_K_rear,
        a1_rear,
        a2_rear,
        args.dt
    )

    alpha_rear = decomposition_degree(
        a1_rear_new,
        a2_rear_new
    )

    ep_rear = np.clip(
        (1.0 - alpha_rear) * EPSILON_V
        + alpha_rear * EPSILON_E,
        0.0,
        1.0,
    )

    save_state_atomic(
        rear_state_file,
        rear_eids,
        a1_rear_new,
        a2_rear_new
    )

    update_serad_loads(
        root,
        rear_eids,
        ep_rear,
        rear_group,
        rear_template,
        "SERad_R_E",
        source_uname=(
            rear_template.get("uname", "")
            if rear_template is not None
            and rear_template.get("uname", "").startswith("Radiation_")
            else None
        ),
    )

    print(
        f"Rear emissivity         : min={ep_rear.min():.6f}, "
        f"max={ep_rear.max():.6f}, mean={ep_rear.mean():.6f}"
    )

    material_uids = build_material_states(root)

    state_by_eid = {
        eid: int(state)
        for eid, state in zip(
            element_ids,
            states_array
        )
    }

    rebuild_state_sets(
        root,
        element_data,
        element_types,
        state_by_eid,
        material_uids
    )

    update_qdec_loads(
        root,
        element_ids,
        qdec,
        qdec_group
    )

    save_state_atomic(
        args.state,
        element_ids,
        a1_new,
        a2_new
    )

    indent_xml(root)
    write_xml_atomic(
        tree,
        args.output
    )

    n_state_elems, n_qdec, n_mats, n_front_serad, n_rear_serad = verify_output(
        args.output,
        element_ids,
        qdec_group,
        expected_front_ids=front_eids,
        expected_rear_ids=rear_eids,
        front_group=front_group,
        rear_group=rear_group,
    )

    print(f"Generated materials : {n_mats}")
    print(f"Generated Q_dec     : {n_qdec}")
    print(f"Generated front Rad : {n_front_serad}")
    print(f"Generated rear Rad  : {n_rear_serad}")
    print(f"State-set elements  : {n_state_elems}")
    print(f"Q_dec group kept    : {qdec_group}")
    print(f"Front Rad group     : {front_group}")
    print(f"Rear Rad group      : {rear_group}")
    print(f"Written             : {args.output}")
    print("XML verification    : OK")
    print("/" * 76)


if __name__ == "__main__":
    main()

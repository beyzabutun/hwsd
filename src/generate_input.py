from gurobipy.gurobipy import multidict


def rack_deployment(num_of_racks, num_of_columns, num_of_groups):
    """
        The info of racks and coordinates and group of racks is generated.

        :param num_of_racks: Number of racks deployed in the WDCN.
        :param num_of_columns: Number of rack columns in the WDCN.
        :param num_of_groups: Number of rack groups in the WDCN.
        :return: (Dict, Dict, Dict)
    """
    rack_info = {}
    r = 0
    c = 0
    order = 0
    for i in range(1, num_of_racks + 1):
        if c >= num_of_columns:
            if num_of_groups == 1:
                order = 0
                c = c - num_of_columns
                r = r + 1
            else:
                order = 1
        if c >= 2 * num_of_columns and order == 1:
            c = c - 2 * num_of_columns
            r = r + 1
            order = 0
        rack_info['R' + str(i)] = [(c, r), order]
        c = c + 1
    racks, rack_coords, rack_order = multidict(rack_info)
    return racks, rack_coords, rack_order


def pm_deployment(num_of_pms, pm_per_rack, num_of_racks, p2u):
    """
    The info of which racks a PM is deployed is generated.

    :param num_of_pms: Number of physical machines in the WDCN.
    :param pm_per_rack: Number of physical machines located in the each rack.
    :param num_of_racks: Number of racks in the WDCN.
    :param p2u: Set of PM to racks.
    :return: Dict
    """
    phy_depl = {}
    ppr = 0
    r = 1
    for pu in p2u:
        phy_depl[pu] = 0
    for p in range(1, num_of_pms + 1):
        if ppr == pm_per_rack:
            ppr = 0
            r = r + 1
        if r > num_of_racks:
            break
        phy_depl['P' + str(p), 'R' + str(r)] = 1
        ppr = ppr + 1
    return phy_depl


def pm_features(num_of_pms):
    """
    List of physical machines and their capacity allowed to be used are generated.
    RAM in GB, cache memory in MB, CPU in core.

    :param num_of_pms: Number of physical machines in the WDCN.
    :return: (Dict, Dict, Dict)
    """
    pm_info = {}
    for p in range(1, num_of_pms+1):
        pm_info['P'+str(p)] = [256, 20, 8]
    phy_machines, allowed_mem_usage, allowed_cache_mem_usage, allowed_cpu_usage = multidict(pm_info)
    return phy_machines, allowed_mem_usage, allowed_cache_mem_usage, allowed_cpu_usage


from generate_input import *
import math
import time
import sys

#PARAMETERS
#Antenna gains Gr,Gt (9 dBi)
Gr = 7.94
Gt = 7.94
Pt = 0.1  #in watt
f = 60.48*math.pow(10,9)  #Hz
c = 3*math.pow(10,8)  #m/s^2
T = 290  #in Kelvin


def calculate_rack_distance(rack, n_of_racks, rack_order, rack_coords, x_dist, y_dist):
    """
        Distances from the given rack to another racks are found.

        :param n_of_racks: Number of racks deployed in the WDCN
        :param rack_order: Holds the group number of each rack in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two racks in x coord.
        :param y_dist: Distance between two racks in y coord.
        :return: List
    """
    dist = [0 for _ in range(n_of_racks)]
    for r in range(n_of_racks):
        dist[r] = calculate_distance(rack, 'R' + str(r + 1), rack_order, rack_coords, x_dist, y_dist)
    return dist


def calculate_distance(rack1, rack2, rack_order, rack_coords, x_dist, y_dist):
    """
        Distance between two racks is calculated.

        :param rack_order: Holds the group number of each rack in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two racks in x coord.
        :param y_dist: Distance between two racks in y coord.
        :return: Float
    """
    x_diff = rack_coords[rack1][0] - rack_coords[rack2][0]
    y_diff = rack_coords[rack1][1] - rack_coords[rack2][1]
    if rack_order[rack1] == rack_order[rack2]:
        d_uv = math.sqrt(math.pow(x_dist * x_diff, 2) + math.pow(y_dist * y_diff, 2))
    else:
        d_uv = math.sqrt(math.pow(abs(x_diff)*x_dist+2, 2) + math.pow(y_dist * y_diff, 2))
    return d_uv


def calculate_p_received(rack1, rack2, rack_order, rack_coords, x_dist, y_dist):
    """
        Received power in scalar and desibel is calculated. (Pr, Pr_dB)

        :param rack_order: Holds the group number of each rack in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two racks in x coord.
        :param y_dist: Distance between two racks in y coord.
        :return: (Float, Float)
    """
    pr = Pt * Gr * Gt * (
        math.pow(c / (4 * 3.14 * 60.48*math.pow(10,9) * calculate_distance(rack1, rack2, rack_order, rack_coords, x_dist, y_dist)), 2)
    )
    return pr, 10 * math.log10(pr)


def calculate_capacity(rack1, rack2, rack_order, rack_coords, x_dist, y_dist):
    """
        Channel capacity between two transceivers located in given two racks is calculated in bps.

        :param rack_order: Holds the group number of the racks in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two racks in x coord.
        :param y_dist: Distance between two racks in y coord.
        :return: Float
    """
    d_uv = calculate_distance(rack1, rack2, rack_order, rack_coords, x_dist, y_dist)
    Pr = Pt * Gr * Gt * (math.pow(c / (4 * 3.14 * 60.48*math.pow(10,9) * d_uv), 2))
    SNR = Pr / (1.38 * math.pow(10, -23) * T * 1.632 * math.pow(10, 9))
    C = 1.632 * math.pow(10, 9) * math.log2(1 + SNR)
    return C


def get_link_capacities(u2v, rack_order, rack_coords, x_dist, y_dist):
    """
        Link capacities between rack pairs are calculated in bps.

        :param u2v: Set of racks to racks
        :param rack_order: Holds the group number of the racks in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two racks in x coord.
        :param y_dist: Distance between two racks in y coord.
        :return: Dict
    """
    l_cap = {}
    for (u, v) in u2v:
        if u != v:
            l_cap[u, v] = calculate_capacity(u, v, rack_order, rack_coords, x_dist, y_dist)
    return l_cap


def get_phy_mach_info(phy_machines, allowed_mem_usage, allowed_cache_mem_usage, allowed_cpu_usage):
    """
        A dictionary holding allowed source usage of each physical machine is created.

        :param allowed_mem_usage: Holds the info. of allowed RAM usage in GB of each physical machine in a dict
        :param allowed_cache_mem_usage: Holds the info. of allowed cache memory usage in MB of each physical machine in
        a dict
        :param allowed_cpu_usage: Holds the info. of allowed CPU usage of each physical machine in a dict
        :return: Dict
    """
    phy_info = {}
    for pm in phy_machines:
        phy_info[pm] = [allowed_mem_usage[pm], allowed_cache_mem_usage[pm], allowed_cpu_usage[pm]]
    return phy_info


def choose_next_rack(current_rack, unvisited_racks, num_of_racks, rack_order, rack_coords, x_dist, y_dist):
    """
        A rack which is closest to the given rack (current rack) is chosen.

        :param unvisited_racks: Holds the list of racks which do not have a wireless link with any other rack
        :param num_of_racks: Number of racks deployed in the WDCN
        :param rack_order: Holds the group number of the racks in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two  racks in x coord.
        :param y_dist: Distance between two racks in y coord.
        :return: Int
    """
    rack_dist = calculate_rack_distance('R'+str(current_rack), num_of_racks, rack_order, rack_coords, x_dist, y_dist)
    max_v = sys.maxsize
    next_rack = None
    for rd in range(len(rack_dist)):
        if rd+1 != current_rack and rd+1 in unvisited_racks and rack_dist[rd] < max_v:
            next_rack = rd+1
            max_v = rack_dist[rd]
    return next_rack


def get_traffic_bw_racks(num_of_racks, rv_info, traffic_info):
    """
        Total traffic demand between rack pairs are calculated and returned as dict.

        :param rv_info: A list holding info. of virtual machines which are deployed in each rack
        {R1: [Vx, Vy, ...], R2: [Vz, ...], ...}
        :param traffic_info: Holds the traffic demands between virtual machine pairs in a dict
        :return: Dict
    """
    traffic_bw_racks = {}
    for r1 in range(1, num_of_racks+1):
        for r2 in range(r1+1, num_of_racks+1):
            traffic_bw_racks['R'+str(r1), 'R'+str(r2)] = 0
            traffic_bw_racks['R'+str(r2), 'R'+str(r1)] = 0
            for vm1 in rv_info[r1]:
                for vm2 in rv_info[r2]:
                    traffic_bw_racks['R'+str(r1), 'R'+str(r2)] = traffic_bw_racks['R'+str(r1), 'R'+str(r2)] + traffic_info[vm1, vm2]  #from r1 to r2
                    traffic_bw_racks['R'+str(r2), 'R'+str(r1)] = traffic_bw_racks['R'+str(r2), 'R'+str(r1)] + traffic_info[vm2, vm1]  #from r2 to r1
    return traffic_bw_racks


def get_virt_deployment(num_of_vms, num_of_pms, vp_info):
    """
        For each virtual machine and physical machine pair, if VM is embedded into PM is returned in a dict.

        :param vp_info: A dict holding the info of in which physical machine that each v. machine is deployed.
        {V1: Px, V2: Py, ...}
        :return: Dict
    """
    vm_info = {}
    for vm in range(1, num_of_vms+1):
        vm_str = 'V'+str(vm)
        if vm_str in vp_info:
            pm_into = vp_info[vm_str]
            for pm in range(1, num_of_pms+1):
                pm_str = 'P'+str(pm)
                if pm_into == pm_str:
                    vm_info[vm_str, pm_str] = 1
                else:
                    vm_info[vm_str, pm_str] = 0

    return vm_info


def get_traffic_bw_2racks(rack1, rack2, rv_info, traffic_info):
    """
        Total traffic demand between two racks is calculated.

        :param rv_info: A dict holding the info of virtual machines which are deployed in each rack
        {R1: [Vx, Vy, ...], R2: [Vz, ...], ...}
        :param traffic_info: Holds the traffic demands between virtual machine pairs
        :return: Float
    """
    traffic = 0
    for vm1 in rv_info[int(rack1[1:])]:
        for vm2 in rv_info[int(rack2[1:])]:
            traffic = traffic + traffic_info[vm1, vm2]  #from rack1 to rack2
    return traffic


def get_demand_info(num_of_vms, rv_info, total_demands_bw_vms):
    """
        Total traffic demand of each VM to the racks (to the VMs located in the rack) is calculated for each rack.
        --traffic_2_r = {R1: {V1: , V2: ... Vm:}, R2: {V1: , V2: ... Vm:}, ... RN: {V1: , V2: ... Vm:}}

        :param num_of_vms: Number of virtual machines to be deployed in the WDCN
        :param rv_info: A dict holding info. of virtual machines which are deployed in each rack
        {R1: [Vx, Vy, ...], R2: [Vz, ...], ...}
        :param total_demands_bw_vms: Holds the total traffic demand between virtual machine pairs in a dict
        :return: Dict
    """
    traffic_2_r = {}
    for rack in rv_info:
        traffic_2_r[rack] = {}
        for vir in range(1, num_of_vms + 1):
            demand = 0
            for vm in rv_info[rack]:
                demand = demand + total_demands_bw_vms[vm][vir - 1]
            traffic_2_r[rack]['V' + str(vir)] = demand
    return traffic_2_r


def hwsd(u2v, rack_order, rack_coords, num_of_vms, num_of_pms, num_of_racks, pm_per_rack, total_demand, traffic_demand,
         total_demands_bw_vms, mem_usage, cache_mem_usage, cpu_usage, phy_depl, phy_machines, allowed_mem_usage,
         allowed_cache_mem_usage, allowed_cpu_usage, x_dist, y_dist, hybrid):
    """
        Heuristic solution is found using Heuristic for Wireless Link and Service Deployment (HWSD) in bps.

        :param u2v: Set of racks to racks
        :param rack_order: Holds the group number of each rack in a dict {R1: #of_group, R2: #of_group,...}
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param num_of_vms: Number of virtual machines to be deployed in the WDCN
        :param num_of_pms: Number of physical machines deployed in the WDCN
        :param num_of_racks: Number of racks deployed in the WDCN
        :param pm_per_rack: Number of physical machines located in each rack
        :param total_demand: Holds total traffic demand of each VM (incoming and outgoing) in a list
        :param traffic_demand: Holds the traffic demand info from a VM to another VM for each pairs in bps in a dict
        :param total_demands_bw_vms: Holds the total traffic demand info between VM pairs (incoming and outgoing) in a dict
        :param mem_usage, cache_mem_usage, cpu_usage: Hold the info. of the source demand of VMs in a dict
        :param phy_depl: Holds the info. of in which rack a physical machine is deployed in a dict
        :param allowed_mem_usage, allowed_cache_mem_usage, allowed_cpu_usage: Holds the info. of allowed source usage
        of each physical machine in a dict
        :param hybrid: Whether the problem is hybrid-WSD (1) or pure-WSD (0)
        :param x_dist: Distance between two racks in x coord. in meters
        :param y_dist: Distance between two racks in y coord. in meters
        :return: Float
    """
    start_time = time.time()
    # Get the capacity of links and physical machines
    lnk_capacities = get_link_capacities(u2v, rack_order, rack_coords, x_dist, y_dist)
    pm_capacities = get_phy_mach_info(phy_machines, allowed_mem_usage, allowed_cache_mem_usage, allowed_cpu_usage)
    total_dem_values = []
    wir_links = {}
    # Get all unvisited racks and virtual machines
    unvisited_services = list(range(1, num_of_vms + 1))
    unvisited_racks = list(range(1, num_of_racks + 1))
    # Set values of variables
    feasible_sol = 0
    current_rack = 1
    # Find the closest rack to establish a wireless link
    next_rack = choose_next_rack(current_rack, unvisited_racks, num_of_racks, rack_order, rack_coords, x_dist, y_dist)
    if next_rack is None:
        print("There is no feasible solution.")
    unvisited_racks.remove(current_rack)
    unvisited_racks.remove(next_rack)
    # wir_links: Holds the wireless link establishment information
    wir_links[current_rack] = next_rack
    wir_links[next_rack] = current_rack
    # Find the service to be deployed
    current_service = total_demand.index(max(total_demand)) + 1
    current_service_str = 'V' + str(current_service)
    total_dem_values.append(current_service_str)
    # rack_traffic: Total traffic of each virtual machine to the vms located in the related rack
    # rack_traffic = {R1: {V1: , V2: ... Vm:}, R2: {V1: , V2: ... Vm:}}
    # vm_traffic = {V1: , V2: ... Vm:}
    rack_traffic = {}
    vm_traffic = {}
    # Hybrid or Pure-WSD
    if hybrid:
        for vir in unvisited_services:
            vm_traffic[vir] = 0
    else:
        for rack in (current_rack, next_rack):
            rack_traffic[rack] = {}
            for vir in unvisited_services:
                rack_traffic[rack][vir] = 0
    # Decrease the capacity of physical machine in which virtual machine deployed
    pm_capacities['P1'][0] = pm_capacities['P1'][0] - mem_usage[current_service_str]
    pm_capacities['P1'][1] = pm_capacities['P1'][1] - cache_mem_usage[current_service_str]
    pm_capacities['P1'][2] = pm_capacities['P1'][2] - cpu_usage[current_service_str]

    print("A wireless link between racks ", current_rack, next_rack)
    print("Deployment of ", current_service_str, "to PM ", 'P1')

    # Hold deployment information
    # {V1: Px, V2: Py...}
    vp_info = {}
    # {R1: [Vx, Vy,...], R2: [Vz,...],...}
    rv_info = {}
    # {V1: Rx, V2: Ry...}
    vr_info = {}
    vp_info[current_service_str] = 'P1'
    rv_info[current_rack] = [current_service_str]
    rv_info[next_rack] = []
    vr_info[current_service_str] = current_rack
    if hybrid:
        for vir in unvisited_services:
            vm_traffic[vir] = vm_traffic[vir] + total_demands_bw_vms[current_service_str][vir - 1]
    else:
        for vir in unvisited_services:
            rack_traffic[current_rack][vir] = rack_traffic[current_rack][vir] + \
                                              total_demands_bw_vms[current_service_str][vir - 1]
    count = 0

    while current_service in unvisited_services:
        # Choose the next service to be deployed
        if hybrid:
            traffic_2_crack = list(vm_traffic.values())
        else:
            traffic_2_crack = list(rack_traffic[current_rack].values())
        max_val = max(traffic_2_crack)
        # Choose next service to be deployed
        next_service = traffic_2_crack.index(max_val) + 1
        next_service_str = 'V' + str(next_service)
        if max_val == 0 and count == 5:
            unvisited_services.remove(current_service)
            print("There are un-deployed services: ", unvisited_services)
            break
        if max_val == 0 and len(unvisited_services) != 0:
            count = count + 1
            temp_rack = current_rack
            current_rack = next_rack
            next_rack = temp_rack
            if hybrid:
                traffic_2_crack = list(vm_traffic.values())
            else:
                traffic_2_crack = list(rack_traffic[current_rack].values())
            max_val = max(traffic_2_crack)
            if max_val == 0 and len(unvisited_services) <= 1:
                unvisited_services.remove(current_service)
                break
            elif max_val == 0 and len(unvisited_services) > 1:
                for uvm in unvisited_services:
                    if uvm != current_service:
                        next_service = uvm
                        next_service_str = 'V' + str(next_service)
            else:
                next_service = traffic_2_crack.index(max_val) + 1
                next_service_str = 'V' + str(next_service)
        # Check the capacity of the link established b/w current and next rack
        # to und. whether it is enough to deploy next_service
        rlink_capa = 0
        llink_capa = 0
        extra_traffic = 0
        for d in rv_info[current_rack]:
            rlink_capa = rlink_capa + traffic_demand[d, next_service_str]
            llink_capa = llink_capa + traffic_demand[next_service_str, d]
        if hybrid:
            for d in rv_info[next_rack]:
                extra_traffic = extra_traffic + traffic_demand[d, next_service_str] + traffic_demand[
                    next_service_str, d]
        if lnk_capacities['R' + str(current_rack), 'R' + str(next_rack)] >= rlink_capa and lnk_capacities[
            'R' + str(next_rack), 'R' + str(current_rack)] >= llink_capa:
            pm_found = 0
            # If capacity is enough, check whether there is a PM with enough capacity
            for p, r in phy_depl.keys():
                if r == 'R' + str(next_rack) and phy_depl[p, r]:
                    if pm_capacities[p][0] >= mem_usage[next_service_str] and \
                            pm_capacities[p][1] >= cache_mem_usage[next_service_str] and \
                            pm_capacities[p][2] >= cpu_usage[next_service_str]:
                        pm_found = 1
                        print("Deployment of ", next_service_str, "to PM", p)
                        rv_info[next_rack].append(next_service_str)
                        vp_info[next_service_str] = p
                        vr_info[next_service_str] = next_rack
                        if hybrid:
                            for vir in unvisited_services:
                                vm_traffic[vir] = vm_traffic[vir] + total_demands_bw_vms[next_service_str][vir - 1]
                            vm_traffic[current_service] = 0
                            vm_traffic[next_service] = 0
                        else:
                            for vir in unvisited_services:
                                rack_traffic[next_rack][vir] = rack_traffic[next_rack][vir] + \
                                                               total_demands_bw_vms[next_service_str][vir - 1]

                            for i in (current_rack, next_rack):
                                for j in (current_service, next_service):
                                    rack_traffic[i][j] = 0
                        unvisited_services.remove(current_service)
                        # Update feasible solution
                        feasible_sol = feasible_sol + rlink_capa + llink_capa + extra_traffic
                        # Decrease the capacities of the related links and PMs
                        lnk_capacities['R' + str(current_rack), 'R' + str(next_rack)] = lnk_capacities['R' + str(
                            current_rack), 'R' + str(next_rack)] - rlink_capa
                        lnk_capacities['R' + str(next_rack), 'R' + str(current_rack)] = lnk_capacities['R' + str(
                            next_rack), 'R' + str(current_rack)] - llink_capa
                        pm_capacities[p][0] = pm_capacities[p][0] - mem_usage[next_service_str]
                        pm_capacities[p][1] = pm_capacities[p][1] - cache_mem_usage[next_service_str]
                        pm_capacities[p][2] = pm_capacities[p][2] - cpu_usage[next_service_str]
                        # Change the rack and VM info
                        temp_rack = current_rack
                        current_rack = next_rack
                        next_rack = temp_rack
                        current_service = next_service
                        break
            if pm_found:
                continue
        # In case of the capacity shortage in wir.links or PMs
        unvisited_services.remove(current_service)
        if len(unvisited_racks) == 0 and len(unvisited_services) != 0:
            print("Couldn't deploy these services: ", unvisited_services)
            break
        # Choose new current and next rack
        current_rack = unvisited_racks[0]
        next_rack = choose_next_rack(current_rack, unvisited_racks, num_of_racks, rack_order, rack_coords, x_dist, y_dist)
        if next_rack == None:
            print("Can't establish a link from ", current_rack, " to another rack.")
            rv_info[current_rack] = []
            break
        wir_links[current_rack] = next_rack
        wir_links[next_rack] = current_rack
        print("A wireless link between racks ", current_rack, next_rack)
        unvisited_racks.remove(next_rack)
        unvisited_racks.remove(current_rack)
        for vm in range(1, num_of_vms + 1):
            if vm not in unvisited_services:
                total_demand[vm - 1] = 0
        current_service = total_demand.index(max(total_demand)) + 1
        current_service_str = 'V' + str(current_service)
        total_dem_values.append(current_service_str)
        current_pm = 'P' + str((current_rack - 1) * (pm_per_rack) + 1)
        # Update deployment values
        vp_info[current_service_str] = current_pm
        rv_info[current_rack] = [current_service_str]
        rv_info[next_rack] = []
        vr_info[current_service_str] = current_rack
        count = 0
        rack_traffic = {}
        vm_traffic = {}
        if hybrid:
            for vir in range(1, num_of_vms + 1):
                vm_traffic[vir] = 0
            for vir in unvisited_services:
                vm_traffic[vir] = vm_traffic[vir] + total_demands_bw_vms[current_service_str][vir - 1]
        else:
            for rack in (current_rack, next_rack):
                rack_traffic[rack] = {}
                for vir in range(1, num_of_vms + 1):
                    rack_traffic[rack][vir] = 0
            for vir in unvisited_services:
                rack_traffic[current_rack][vir] = rack_traffic[current_rack][vir] + \
                                                  total_demands_bw_vms[current_service_str][vir - 1]

        # Decrease the capacity of physical machine into which virtual machine deployed
        pm_capacities[current_pm][0] = pm_capacities[current_pm][0] - mem_usage[current_service_str]
        pm_capacities[current_pm][1] = pm_capacities[current_pm][1] - cache_mem_usage[current_service_str]
        pm_capacities[current_pm][2] = pm_capacities[current_pm][2] - cpu_usage[current_service_str]
        print("Deployment of ", current_service_str, "to PM ", current_pm)

    total_time = time.time() - start_time
    print('Found feasible solution, value %g %.4f seconds' % (feasible_sol, total_time))
    # Calculate total demand of each virtual machine to each rack
    traffic_2_r = get_demand_info(num_of_vms, rv_info, total_demands_bw_vms)
    # Calculate all traffic demand between rack pairs
    traffic_bw_racks = get_traffic_bw_racks(num_of_racks, rv_info, traffic_demand)
    # Set 1 if v is in p, otherwise 0. To use in optimization problems
    vm_info = get_virt_deployment(num_of_vms, num_of_pms, vp_info)
    return unvisited_services, unvisited_racks, traffic_2_r, traffic_bw_racks, rv_info, vr_info, \
           vp_info, vm_info, pm_capacities, lnk_capacities, wir_links, feasible_sol


def i_hwsd(num_of_vms, num_of_racks, rack_order, rack_coords, x_dist, y_dist, phy_depl,
           unvisited_vms, unvisited_racks, traffic_2_r, traffic_demand, vr_info, rv_info, vp_info, vm_info,
           pm_capacity, link_capacity, link_info, mem_usage, cache_mem_usage, cpu_usage, feasible_sol,
           total_demands_bw_vms, hybrid):

    """
        Heuristic solution found by HWSD is improved using Improved Heuristic for Wireless Link and Service Deployment (I-HWSD)
        in bps.

        :param rack_order: Holds the group number of each rack in a dict
        :param rack_coords: Holds the coord. info. of the racks in a dict
        :param x_dist: Distance between two racks in x coord. in meters
        :param y_dist: Distance between two racks in y coord. in meters
        :param phy_depl: Holds the info. of in which rack a physical machine is deployed ln a dict
        :param traffic_demand: Holds the traffic demand info from a VM to another VM for each pairs in bps in a dict
        :param total_demands_bw_vms: Holds the total traffic demand info between VM pairs (incoming and outgoing) in a dict
        :param mem_usage, cache_mem_usage, cpu_usage: Hold the info. of the source demand of VMs in a dict
        :param hybrid: Whether the problem is hybrid-WSD (1) or pure-WSD (0)
        ## THE FOLLOWINGS ARE RETURNED FROM HWSD
        :param unvisited_vms: List of un-deployed VMs in hwsd
        :param unvisited_racks: List of racks which have no wireless link established in hwsd
        :param traffic_2_r: Contains the total traffic demand of each virtual machine to each rack in a dict
        ={R1: {V1: , V2:, ... Vm}, ...}
        :param vr_info: Contains the VM and rack pairs where the VM is deployed into the rack in a dict
        :param rv_info: A dict holding virtual machines which are deployed in each rack
        {R1: [Vx, Vy, ...], R2: [Vz, ...], ...}
        :param vp_info: A dict holding the info of in which physical machine each v. machine is deployed.
        {V1: Px, V2: Py, ...}
        :param vm_info: Contains the info of in which PM a VM is deployed in a dict
        {(Vx, Py): 1 or 0...}
        :param pm_capacity: Holds the remaining capacity (from hwsd) of each PM in terms of RAM, cache memory, and CPU
        in a dict
        :param link_capacity: Shows the remaining capacity (from hwsd) of each wireless link established in a dict
        :param link_info: Shows the rack pairs between which a wireless link is established in hwsd
        :param feasible_sol: Solution obtained by HWSD.
        :return: Float
    """
    start_time = time.time()
    epsilon = 5.
    total_gain = 5.1
    # Continue while total traffic demand gain > epsilon value
    while total_gain > epsilon or len(unvisited_vms) != 0:
        total_gain = 0.
        # Repeat for all VMs
        for vm in range(1, num_of_vms + 1):
            vm_str = 'V' + str(vm)
            unvisited = 1
            # Check whether the VM is deployed in hwsd algorithm
            # Save info. of the rack in which the related VM deployed
            if vm not in unvisited_vms:
                current_rack = vr_info[vm_str]
                try:
                    neighbor_rack = link_info[current_rack]
                except:
                    continue
                current_pm = vp_info[vm_str]
                unvisited = 0
                # Calculate the traffic demand contribution of the VM to the total traffic demand satisfied
                if hybrid:
                    contr = traffic_2_r[neighbor_rack][vm_str] + traffic_2_r[current_rack][vm_str]
                else:
                    contr = traffic_2_r[neighbor_rack][vm_str]
            else:
                contr = -1
            deployment_info = [0 for _ in range(0, 7)]
            # For all other racks,
            # Calculate the traffic demand contribution of the related VM
            # in case of deploying it in different racks
            for rack in rv_info:
                try:
                    rack_into = link_info[rack]
                except:
                    if len(unvisited_racks) > 1:
                        rack_into = choose_next_rack(rack, unvisited_racks, num_of_racks, rack_order, rack_coords,
                                                     x_dist, y_dist)
                        if rack_into is None:
                            break
                        link_info[rack] = rack_into
                        link_info[rack_into] = rack
                        unvisited_racks.remove(rack)
                        unvisited_racks.remove(rack_into)
                    else:
                        break
                if hybrid:
                    demand = traffic_2_r[rack][vm_str] + traffic_2_r[rack_into][vm_str]
                else:
                    demand = traffic_2_r[rack][vm_str]
                # Consider re-deployment decision
                if demand > contr and demand - contr > deployment_info[4]:
                    # Check capacity
                    # If it is enough -> deploy -> calculate gain -> calculate new solution
                    # If not, continue
                    rlink_capa = 0
                    llink_capa = 0
                    for d in rv_info[rack]:
                        if d != vm_str:
                            rlink_capa = rlink_capa + traffic_demand[vm_str, d]
                            llink_capa = llink_capa + traffic_demand[d, vm_str]
                    for p, r in phy_depl.keys():
                        if r == 'R' + str(rack_into) and phy_depl[p, r]:
                            if pm_capacity[p][0] >= mem_usage[vm_str] and \
                                    pm_capacity[p][1] >= cache_mem_usage[vm_str] and \
                                    pm_capacity[p][2] >= cpu_usage[vm_str]:
                                deployment_info[0] = rack_into
                                deployment_info[1] = rack
                                deployment_info[2] = p
                                deployment_info[3] = demand
                                if not unvisited:
                                    deployment_info[4] = demand - contr  # gain
                                if unvisited:
                                    contr = 0
                                    deployment_info[4] = demand - contr  # gain
                                deployment_info[5] = rlink_capa
                                deployment_info[6] = llink_capa
                                break
            # Update after re-deployment decision
            if deployment_info[0] != 0:
                vm_change = vm_change + 1
                total_gain = total_gain + deployment_info[4]
                feasible_sol = feasible_sol + deployment_info[4]
                rack_into = deployment_info[0]
                rack = deployment_info[1]
                pm = deployment_info[2]
                rlink_capa = deployment_info[5]
                llink_capa = deployment_info[6]
                if unvisited:
                    print("old rack:", "It was an unvisited service")
                    unvisited_vms.remove(vm)
                else:
                    print("old rack:", current_rack)
                print("new rack:", rack_into)
                print("gain:", deployment_info[4])
                print("SOL: ", feasible_sol)
                # Decrease the capacities of the related links and PMs
                link_capacity['R' + str(rack_into), 'R' + str(rack)] = link_capacity['R' + str(rack_into), 'R' + str(
                    rack)] - rlink_capa
                link_capacity['R' + str(rack), 'R' + str(rack_into)] = link_capacity['R' + str(rack), 'R' + str(
                    rack_into)] - llink_capa
                pm_capacity[pm][0] = pm_capacity[pm][0] - mem_usage[vm_str]
                pm_capacity[pm][1] = pm_capacity[pm][1] - cache_mem_usage[vm_str]
                pm_capacity[pm][2] = pm_capacity[pm][2] - cpu_usage[vm_str]
                rv_info[rack_into].append(vm_str)
                if not unvisited:
                    rv_info[current_rack].remove(vm_str)
                    vm_info[vm_str, current_pm] = 0
                vr_info[vm_str] = rack_into
                vm_info[vm_str, pm] = 1
                vp_info[vm_str] = pm
                # Get new traffic_2_r after changes : {R1: {V1: ... Vm:}, ..., Rn:{V1:, ..., V2:}}
                traffic_2_r = get_demand_info(num_of_vms, rv_info, total_demands_bw_vms)
                ## UPDATE OLD DEPLOYMENT VALUES AFTER RE-DEPL. DECISION ##
                if not unvisited:
                    rlink_capa = 0
                    llink_capa = 0
                    for d in rv_info[neighbor_rack]:
                        if d != vm_str:
                            rlink_capa = rlink_capa + traffic_demand[vm_str, d]
                            llink_capa = llink_capa + traffic_demand[d, vm_str]
                    link_capacity['R' + str(current_rack), 'R' + str(neighbor_rack)] = link_capacity['R' + str(
                        current_rack), 'R' + str(neighbor_rack)] + rlink_capa
                    link_capacity['R' + str(neighbor_rack), 'R' + str(current_rack)] = link_capacity['R' + str(
                        neighbor_rack), 'R' + str(current_rack)] + llink_capa
                    pm_capacity[current_pm][0] = pm_capacity[current_pm][0] + mem_usage[vm_str]
                    pm_capacity[current_pm][1] = pm_capacity[current_pm][1] + cache_mem_usage[vm_str]
                    pm_capacity[current_pm][2] = pm_capacity[current_pm][2] + cpu_usage[vm_str]
      total_time = time.time() - start_time
      return feasible_sol, total_time

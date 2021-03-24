import csv

from scipy.spatial import distance

import Parameter as para
from Network_Method import uniform_com_func, to_string, count_package_function, partition_function, get_D_max, get_ECR_max, get_CN_max


class Network:
    def __init__(self, list_node=None, mc_list = None, target=None):
        self.node = list_node
        self.set_neighbor()
        self.set_level()
        self.mc_list = mc_list
        self.target = target
        self.began = False
        self.D_avg = []
        self.request_id = []
        self.D_max = 0
        self.ECR_max = 0

    def set_angle(self):
        for node in self.node:
            node.set_angle()

    def set_neighbor(self):
        for node in self.node:
            for other in self.node:
                if other.id != node.id and distance.euclidean(node.location, other.location) <= node.com_ran:
                    node.neighbor.append(other.id)

    def set_level(self):
        queue = []
        for node in self.node:
            if distance.euclidean(node.location, para.base) < node.com_ran:
                node.level = 1
                queue.append(node.id)
        while queue:
            for neighbor_id in self.node[queue[0]].neighbor:
                if not self.node[neighbor_id].level:
                    self.node[neighbor_id].level = self.node[queue[0]].level + 1
                    queue.append(neighbor_id)
            queue.pop(0)

    def communicate(self, func=uniform_com_func):
        return func(self)

    # def partition(self, partition_func=partition_function):
    #     partition_func(self)
    #     self.partitioned = True


    def run_per_second(self, t, optimizer):
        state = self.communicate()
        self.request_id = []
        if self.began:   
            if optimizer:
                all_mc_deactivated = True
                for mc in self.mc_list:
                    if mc.is_active == True:
                        all_mc_deactivated = False
                if all_mc_deactivated:
                    optimizer.get_charging_route(network=self)
                for mc in self.mc_list:
                    mc.run(network=self, time_stem=t, optimizer=optimizer)
            return state

    def simulate_lifetime(self, optimizer, file_name="log/energy_log.csv", D_max_func=get_D_max, ECR_max_func=get_ECR_max):
        energy_log = open(file_name, "w")
        writer = csv.DictWriter(energy_log, fieldnames=["time", "mc energy", "min energy"])
        writer.writeheader()
        t = 0
        while self.node[self.find_min_node()].energy >= 0 and t <= 2000000:
            t = t + 1
            if (t-1) % 100 == 0:
                print("time = ", t, ", lowest enegy node: ", self.node[self.find_min_node()].energy, "at", self.node[self.find_min_node()].location)
                for mc in self.mc_list:
                    print("\tMC#{} at{} is {}".format(mc.id, mc.current, mc.get_status()))
            state = self.run_per_second(t, optimizer)
            if t == 200:
                for index, node in enumerate(self.node):
                    if (t - node.check_point[-1]["time"]) > 50:
                        node.set_check_point(t)
                self.began = True 
        
        print(t, self.node[self.find_min_node()].energy)
        # for mc in self.mc_list:
        #     print("\tMC#{} at{}".format(mc.id, mc.current))
        #     writer.writerow({"time": t, "mc energy": mc.energy, "min energy": self.node[self.find_min_node()].energy})
        energy_log.close()
        return t

    def simulate_max_time(self, optimizer, max_time=10000, file_name="log/information_log.csv"):
        information_log = open(file_name, "w")
        writer = csv.DictWriter(information_log, fieldnames=["time", "nb dead", "nb package"])
        writer.writeheader()
        nb_dead = 0
        nb_package = len(self.target)
        t = 0
        while t <= max_time:
            t += 1
            if (t-1)%100 == 0:
                print(t, self.mc.current, self.node[self.find_min_node()].energy)
            state = self.run_per_second(t, optimizer)
            current_dead = self.count_dead_node()
            current_package = self.count_package()
            if current_dead != nb_dead or current_package != nb_package:
                nb_dead = current_dead
                nb_package = current_package
                writer.writerow({"time": t, "nb dead": nb_dead, "nb package": nb_package})
        print(t, self.mc.current, self.node[self.find_min_node()].energy)
        information_log.close()
        return t

    def simulate(self, optimizer, max_time=None, file_name="log/energy_log.csv"):
        if max_time:
            life_time = self.simulate_max_time(optimizer=optimizer, max_time=max_time, file_name=file_name)
        else:
            life_time = self.simulate_lifetime(optimizer=optimizer, file_name=file_name)
        return life_time

    def print_net(self, func=to_string):
        func(self)

    def find_min_node(self):
        min_energy = 10 ** 10
        min_id = -1
        for node in self.node:
            if node.energy < min_energy:
                min_energy = node.energy
                min_id = node.id
        return min_id

    def count_dead_node(self):
        count = 0
        for node in self.node:
            if node.energy < 0:
                count += 1
        return count

    def count_package(self, count_func=count_package_function):
        count = count_func(self)
        return count

    def request(self, mc, t):
        region = mc.id
        E_max = self.node[0].energy_max
        N_coverage = (mc.energy - self.D_avg[region]*mc.e_move/mc.velocity) / (self.D_avg[region]*mc.e_move/mc.velocity + E_max)
        T_coverage = (N_coverage + 1)*(self.D_avg[region]/mc.velocity) + N_coverage*E_max/0.04
        for node in self.partitioned_node[region]:
            if (1000 - node.check_point[-1]["time"]) > 50:
                node.set_check_point(t)
            C_thresh = node.avg_energy*T_coverage*len(self.partitioned_node[region])/N_coverage
            if node.energy < C_thresh:
                node.request(mc, t)
                self.request_id.append(node.id)
    def get_max_D(self, D_max_func=get_D_max):
        self.D_max = D_max_func(self)
        return self.D_max

    def get_max_ECR(self, ECR_max_func=get_ECR_max):
        self.ECR_max = ECR_max_func(self)
        return self.ECR_max
    
    def get_max_CN(self, CN_max_func=get_CN_max):
        return CN_max_func(self)
    

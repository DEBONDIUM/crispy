#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 17:41:12 2025

@author: lbremaud
"""

# =============================================================================
# LIB
# =============================================================================
from typing import List, Union
import open3d as o3d
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import os
import glob
import networkx as nx
import random
import sys

# =============================================================================
# CLASS Viewer3D
# =============================================================================
class Viewer3D:
    def __init__(self, crispy_instance):
        self.crispy = crispy_instance
        self.iteration = 0
        self.vis = o3d.visualization.VisualizerWithKeyCallback()
        self.cloud = o3d.geometry.PointCloud()

    def _update_geometry(self):
        fragments = self.crispy.fragments[self.iteration]
        points = []
        colors = []
        for fragment in fragments:
            if fragment.ID != -1:
                pts = fragment.nodes[:, :3]
                color = np.array(self.crispy._random_color_ID(fragment.ID)).reshape(1, 3)
                points.append(pts)
                colors.append(np.tile(color, (pts.shape[0], 1)))
        if points:
            pts = np.vstack(points)
            cols = np.vstack(colors)
            self.cloud.points = o3d.utility.Vector3dVector(pts)
            self.cloud.colors = o3d.utility.Vector3dVector(cols)

        self._print_progress_bar(self.iteration, len(self.crispy.fragments))

    def _next(self, vis):
        if self.iteration < len(self.crispy.fragments) - 1:
            self.iteration += 1
            self._update_geometry()
            vis.update_geometry(self.cloud)
            vis.poll_events()
            vis.update_renderer()
        return False

    def _prev(self, vis):
        if self.iteration > 0:
            self.iteration -= 1
            self._update_geometry()
            vis.update_geometry(self.cloud)
            vis.poll_events()
            vis.update_renderer()
        return False

    def _random_color_ID(self, ID):
        np.random.seed(ID)
        return np.random.rand(3)

    def _print_progress_bar(self, iteration, total, length=40):
        percent = iteration / (total - 1) if total > 1 else 1
        filled = int(length * percent)
        bar = '█' * filled + '-' * (length - filled)
        text = f"Progress: |{bar}| {iteration + 1}/{total}"
    
        sys.stdout.write('\r' + text.ljust(80))
        sys.stdout.flush()

    def run(self):
        self.vis.create_window(window_name="Crispy 3D Viewer", width=1024, height=768)
        self._update_geometry()
        self.vis.add_geometry(self.cloud)
    
        self.vis.register_key_callback(263, self._prev)
        self.vis.register_key_callback(262, self._next)
    
        self.vis.run()
        self.vis.destroy_window()
        print()

# =============================================================================
# CLASS FragmentDetector
# =============================================================================
class FragmentDetector:
    def __init__(self, source:Union[str, List[str]], trash_treshold:int = 10):
        self._trash_treshold = trash_treshold

        if isinstance(source, str):
            self._directory = source
            self._agdd_files = self._load_agdd()
        elif isinstance(source, list):
            self._directory = os.path.commonpath(source)
            self._agdd_files = sorted(source)
        else:
            raise ValueError("source must be a directory path or a list of .agdd files")
        
        self._fragments = []
        self._iterations_nb = len(self._agdd_files)
        self._seed = random.randint(0, 2**32 - 1)
        self._2d_plot_lim = []
    
    def build_fragments(self):
        for iteration, agdd_file in enumerate(self._agdd_files):
            new_fragments = Fragment.detect_from_agdd(self, agdd_file, iteration, self._trash_treshold)
            self._fragments.append(new_fragments)

    def _load_agdd(self):
        agdd_files = glob.glob(os.path.join(self._directory, "*.agdd"))
        agdd_files.sort()
        return agdd_files
    
    def get_fragments_by_ID(self, ID:int):
        list_ = []
        for fragment_list in self._fragments:
            for fragment in fragment_list:
                if fragment.ID == ID:
                    list_.append(fragment)
        return list_
    
    def _get_fragments_by_iteration_safe(self, iteration:int):
        if 0 <= iteration < len(self._fragments):
            return self._fragments[iteration]
        else:
            return []

    def get_fragments_by_iteration(self, iteration:int):
        if iteration < 0:
            raise ValueError("Iteration must be non-negative.")
        if iteration >= len(self._fragments):
            raise IndexError(f"Iteration {iteration} is out of range (max is {len(self._fragments) - 1}).")
        return self._fragments[iteration]
    
    def save(self):
        if not self.fragments:
            raise ValueError("No fragments built!")
            
        os.makedirs(os.path.join(self._directory, 'txt'), exist_ok=True)
            
        for fragment_list in self._fragments:
            for fragment in fragment_list:
                ID = fragment.ID
                iteration = fragment.iteration
                area = fragment.area
                parentIDs = fragment.parentIDs
                parentIDs = '\t'.join(str(pid) for pid in fragment.parentIDs)
                nodeIDs = fragment.nodeIDs
                nodes = fragment.nodes
                bondIDs = fragment.bondIDs
                bonds = fragment.bonds

                
                filename = os.path.join(self._directory, 'txt', f"fragment_ID{ID}_iteration{iteration}.txt")
                with open(filename, 'w') as file:
                    file.write(
                        f"# ID\n{ID}\n"
                        f"# Iteration\n{iteration}\n"
                        f"# Area [%]\n{area}\n"
                         "# Parent IDs\n"
                    )
                    file.write(parentIDs + "\n")
                
                    # Write node IDs and their coordinates
                    file.write("# NodeID / X / Y / Z / R\n")
                    for node_id, (x, y, z, r) in zip(nodeIDs, nodes):
                        file.write(f"{node_id}\t{x:.6f}\t{y:.6f}\t{z:.6f}\t{r:.6f}\n")
                
                    # Write bond IDs and the nodes they connect
                    file.write("# BondID / Node1 / Node2\n")
                    for bond_id, (n1, n2) in zip(bondIDs, bonds):
                        file.write(f"{bond_id}\t{n1}\t{n2}\n")
        print(f"Fragments data have been saved in {os.path.join(self._directory, 'txt')}.")
    
    @property
    def directory(self):
        return self._directory
    @property
    def trash_treshold(self):
        return self._trash_treshold
    @property
    def agdd_files(self):
        return self._agdd_files
    @property
    def fragments(self):
        return self._fragments
    @property
    def iterations_nb(self):
        return self._iterations_nb
    
    def _random_color_ID(self, ID:int):
        np.random.seed(self._seed + ID)
        r, g, b = np.random.rand(3)
        return (r, g, b)
    
    def plot2D(self,
               iteration:int,
               show_trash:bool    = False,
               show_bonds:bool    = False,
               auto_close:bool    = False,
               dot_size:int       = 4,
               rigid_lim:bool     = True,
               border_coef:float  = 1.5,
               save:bool          = True,
               save_format:str    = "png"):
        
        if not self.fragments:
            raise ValueError("No fragments built!")
        
        fragments = self.get_fragments_by_iteration(iteration)

        fig, ax = plt.subplots(figsize = (6, 6))
        for fragment in fragments:
            nodes = fragment.nodes
            ID = fragment.ID
            nodeIDs = fragment.nodeIDs
            bonds = fragment.bonds
            
            if ID != -1:
                ax.scatter(nodes[:, 0], nodes[:, 1], s=dot_size, c=[self._random_color_ID(ID)])
                centroid = np.mean(nodes, axis=0)
                ax.text(centroid[0], centroid[1], f"{ID}", color='k', fontsize=10, weight='bold', clip_on=True)
                
                if show_bonds:
                    for bond in bonds:
                        x1 = nodes[np.where(bond[0] == nodeIDs)[0], 0]
                        y1 = nodes[np.where(bond[0] == nodeIDs)[0], 1]
                        x2 = nodes[np.where(bond[1] == nodeIDs)[0], 0]
                        y2 = nodes[np.where(bond[1] == nodeIDs)[0], 1]
                        ax.plot([x1, x2], [y1, y2], '-', c=self._random_color_ID(ID))
            elif show_trash:
                ax.scatter(nodes[:, 0], nodes[:, 1], s=dot_size, c='k')
        
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(labelbottom=False, labelleft=False)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.axis('equal')
        if rigid_lim and iteration == 0:
            self._2d_plot_lim = [(ax.get_xlim()[0] * border_coef, ax.get_xlim()[1] * border_coef),
                                 (ax.get_ylim()[0] * border_coef, ax.get_ylim()[1] * border_coef)]
        ax.set_xlim(self._2d_plot_lim[0])
        ax.set_ylim(self._2d_plot_lim[1])
        plt.tight_layout()
        
        if save:
            os.makedirs(os.path.join(self._directory, 'img'), exist_ok=True)
            filename = os.path.join(self._directory, 'img', f"plot2D_iteration{iteration}.{save_format}")
            plt.savefig(filename, format=save_format)
        if auto_close:
            plt.close(fig)
    
    def stackplot(self,
                  xgrid:bool         = True,
                  interline:bool     = True,
                  auto_close:bool    = False,
                  save:bool          = True,
                  save_format:str    = "png"):
        
        if not self.fragments:
            raise ValueError("No fragments built!")
                
        IDs = []
        for fragment_list in self._fragments:
            for fragment in fragment_list:
                IDs.append(fragment.ID)
        IDs = np.unique(IDs)
        
        colors = ['k'] + [self._random_color_ID(ID) for ID in IDs if ID != -1]
        iterations = np.arange(self.iterations_nb, dtype=int)
        
        data = np.zeros((len(IDs), len(iterations)))
        for i, ID in enumerate(IDs):
            fragments = self.get_fragments_by_ID(ID)
            for fragment in fragments:
                data[i, fragment.iteration] = fragment.area
        cumulative_data = np.cumsum(data, axis=0)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.stackplot(iterations, data, labels=IDs, colors=colors, alpha=0.8)
        if xgrid:
            for iteration in iterations: ax.axvline(x=iteration, color='k', linestyle='--', linewidth=0.5)
        if interline:
            for line in cumulative_data[:-1]: ax.plot(iterations, line, color="black", linewidth=1)
        
        ax.legend(loc='center left', framealpha=1., bbox_to_anchor=(1., 0.5))
        ax.set_xlim(0, np.max(iterations) - 1)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Iterations", fontsize=10)
        ax.set_ylabel("Area / %", fontsize=10)
        plt.tight_layout()
        
        if save:
            os.makedirs(os.path.join(self._directory, 'img'), exist_ok=True)
            filename = os.path.join(self._directory, 'img', f"stackplot.{save_format}")
            plt.savefig(filename, format=save_format)
        if auto_close:
            plt.close(fig)
    
    def graphplot(self,
                  auto_close:bool    = False,
                  save:bool          = True,
                  save_format:str    = "png"):
        
        if not self.fragments:
            raise ValueError("No fragments built!")

        IDs = []
        for fragment_list in self._fragments:
            for fragment in fragment_list:
                IDs.append(fragment.ID)
        IDs = np.unique(IDs)

        tree = []
        for ID in IDs:
            if ID != -1:
                for fragment in  self.get_fragments_by_ID(ID):
                    tree.append([ID, fragment.parentIDs, fragment.iteration, fragment.area])
    
        G = nx.DiGraph()

        areas = {}
        for fragment in tree:
            ID, parentIDs, iteration, area = fragment
            areas[ID] = area
            for parentID in parentIDs:
                if parentID != ID:
                    G.add_edge(parentID, ID)
        
        pos = nx.nx_agraph.graphviz_layout(G, prog="twopi")
            
        min_size = 500
        sizes = [max(areas[node] * min_size, min_size) for node in G.nodes]
        
        # colors
        colors = [self._random_color_ID(node) for node in G.nodes]
    
        fig = plt.figure(figsize=(12, 8))
        nx.draw(G, pos,
                with_labels=True,
                node_size=sizes,
                node_color=colors,
                font_size=10,
                font_weight="bold",
                edge_color="k",
                edgecolors="k")
        
        if save:
            os.makedirs(os.path.join(self._directory, 'img'), exist_ok=True)
            filename = os.path.join(self._directory, 'img', f"graphplot.{save_format}")
            plt.savefig(filename, format=save_format)
        if auto_close:
            plt.close(fig)
    
    def viewer3D(self):
        if not self.fragments:
            raise ValueError("No fragments built!")
            
        viewer = Viewer3D(self)
        viewer.run()

# =============================================================================
# CLASS Fragment
# =============================================================================
class Fragment:
    def __init__(self, ID:int, nodeIDs:list, nodes:list, bondIDs:list, bonds:list, iteration:int, parentIDs:list, area:float):
        self._ID = ID
        self._nodeIDs = nodeIDs
        self._nodes = nodes
        self._bondIDs = bondIDs
        self._bonds = bonds
        self._iteration = iteration
        self._parentIDs = parentIDs
        self._area = area
    
    @property
    def ID(self):
        return self._ID
    @property
    def nodeIDs(self):
        return self._nodeIDs
    @property
    def nodes(self):
        return self._nodes
    @property
    def bondIDs(self):
        return self._bondIDs
    @property
    def bonds(self):
        return self._bonds
    @property
    def iteration(self):
        return self._iteration
    @property
    def parentIDs(self):
        return self._parentIDs
    @property
    def area(self):
        return self._area
    
    @classmethod
    def detect_from_agdd(cls, detector:'FragmentDetector', agdd_file:str, iteration:int, trash_treshold:int):
        nodes_nb, bonds_nb, nodes_all, bonds_all = cls._read_agdd(agdd_file)
        connectivity = cls._compute_connectivity(nodes_nb, bonds_all)
        fragments = cls._identify_fragments(detector, connectivity, nodes_all, bonds_all, iteration, trash_treshold)
        return fragments

    @staticmethod
    def _read_agdd(agdd_file:str):
        file = open(agdd_file, "r")
        lines = file.readlines()
        file.close()
        
        nodes_nb = int(lines[0])
        bonds_nb = int(lines[nodes_nb + 1])
        
        nodes_all = np.array([list(map(float, line.split())) for line in lines[1:nodes_nb + 1]])
        bonds_all = np.array([list(map(int, line.split())) for line in lines[nodes_nb + 2:nodes_nb + 2 + bonds_nb]])
        return nodes_nb, bonds_nb, nodes_all, bonds_all
    
    @staticmethod
    def _compute_connectivity(nodes_nb:int, bonds_all:list):
        connectivity = lil_matrix((nodes_nb, nodes_nb), dtype=int)
        for i, j in bonds_all:
            connectivity[i, j] = 1
            connectivity[j, i] = 1
        return connectivity.tocsr()
    
    def _identify_fragments(detector:'FragmentDetector', connectivity:csr_matrix, nodes_all:list, bonds_all:list, iteration:int, trash_treshold:int):
        # initialize fragment list
        fragments = []
        fragmentIDs = []
        
        # get previous fragments
        fragments_previous = detector._get_fragments_by_iteration_safe(iteration - 1)
        
        # compute area_all
        area_all = np.pi * np.sum(nodes_all[:, 3] ** 2.)
        
        # compute fragments
        fragments_nb, nodes_label_all = connected_components(connectivity, directed=False)
        
        # fragments labels and nodes : divide by trash_treshold size
        labels_and_nodeIDs_all = np.empty((fragments_nb, 2), dtype=object)
        for i in range(fragments_nb):
            labels_and_nodeIDs_all[i, 0] = i
            labels_and_nodeIDs_all[i, 1] = np.where(nodes_label_all == i)[0]
            
        labels_and_nodeIDs_trash = np.array([f for f in labels_and_nodeIDs_all if len(f[1]) <= trash_treshold])
        labels_and_nodeIDs_fragments = np.array([f for f in labels_and_nodeIDs_all if len(f[1]) > trash_treshold])
        
        # compute trash
        if len(labels_and_nodeIDs_trash) > 0:
            # get nodes
            nodeIDs_trash = np.concatenate(labels_and_nodeIDs_trash[:, 1])
            nodes_trash = nodes_all[nodeIDs_trash]
            
            # get bonds
            mask = np.isin(bonds_all[:, 0], nodeIDs_trash) | np.isin(bonds_all[:, 1], nodeIDs_trash)
            bondIDs_trash = np.where(mask)[0]
            bonds_trash = bonds_all[bondIDs_trash]
            
            # get parentIDs
            parentIDs_trash = []
            for fragment in fragments_previous:
                if fragment.ID != -1:
                    mask = np.isin(nodeIDs_trash, fragment.nodeIDs)
                    if np.any(mask):
                        parentIDs_trash.append(fragment.ID)
            
            # compute area partition
            area_trash = 100. * (np.pi * np.sum(nodes_trash[:, 3] ** 2.)) / area_all

            # create fragment
            fragments.append(Fragment(-1, nodeIDs_trash, nodes_trash, bondIDs_trash, bonds_trash, iteration, sorted(parentIDs_trash), area_trash))

        # compute fragments
        for i, labels_and_nodeIDs in enumerate(labels_and_nodeIDs_fragments):
            # get nodes
            nodeIDs_fragment = labels_and_nodeIDs_fragments[i, 1]
            nodes_fragment = nodes_all[nodeIDs_fragment]
            
            # get bonds
            mask = np.isin(bonds_all[:, 0], nodeIDs_fragment) | np.isin(bonds_all[:, 1], nodeIDs_fragment)
            bondIDs_fragment = np.where(mask)[0]
            bonds_fragment = bonds_all[bondIDs_fragment]
            
            # get parents
            parentIDs_fragment = []
            if iteration > 0:
                for fragment in fragments_previous:
                    if nodeIDs_fragment[0] in fragment.nodeIDs:
                        parentIDs_fragment = [fragment.ID]
                        continue

            # get ID
            ID_fragment = nodeIDs_fragment[0]
            if iteration > 0:
                if parentIDs_fragment[0] not in fragmentIDs:
                    ID_fragment = parentIDs_fragment[0]
                else:
                    ID_fragment = nodeIDs_fragment[0]
            
            # compute area partition
            area_fragment = 100. * (np.pi * np.sum(nodes_fragment[:, 3] ** 2.)) / area_all
            
            # create fragments list
            fragments.append(Fragment(ID_fragment, nodeIDs_fragment, nodes_fragment, bondIDs_fragment, bonds_fragment, iteration, sorted(parentIDs_fragment), area_fragment))
            fragmentIDs.append(ID_fragment)
        
        return fragments
    
# =============================================================================
# MAIN
# =============================================================================
# sim = Crispy(["../examples/agdd/domain-0000000004000.agdd", "../examples/agdd/domain-0000000005000.agdd"])
# sim = Crispy("../examples/agdd/")
# sim.build_fragments()
# # sim.plot3D()
# sim.save()

# # save= True
# # for i in range(sim.iterations_nb):
# #     sim.plot2D(iteration=i, save = save)

# # sim.stackplot()
# # sim.graphplot()

# sim.viewer3D()





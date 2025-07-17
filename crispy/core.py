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
        """
        3D viewer of the fragments over all iterations.
    
        Args:
            crispy_instance (FragmentDetector): An instance containing fragments and related data.
        """
        self.crispy = crispy_instance
        self.iteration = 0
        self.vis = o3d.visualization.VisualizerWithKeyCallback()
        self.cloud = o3d.geometry.PointCloud()

# === Core methods ===
    def run(self) -> None:
        """
        Run the 3D viewer.
        """
        self.vis.create_window(window_name="Crispy 3D Viewer", width=1024, height=768)
        self._update_geometry()
        self.vis.add_geometry(self.cloud)
    
        self.vis.register_key_callback(263, self._prev)
        self.vis.register_key_callback(262, self._next)
    
        self.vis.run()
        self.vis.destroy_window()
        print()

# === Internal helpers ===
    def _update_geometry(self) -> None:
        """
        Update the 3D viewer window.
        """
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

    def _next(self,
        vis:o3d.visualization.Visualizer
        ) -> bool:
        """
        Callback to go to the next iteration in the 3D viewer.
    
        Args:
            vis (open3d.visualization.Visualizer): The Open3D visualizer instance triggering the callback.
    
        Returns:
            bool: False to indicate that the visualizer should continue running.
        """
        if self.iteration < len(self.crispy.fragments) - 1:
            self.iteration += 1
            self._update_geometry()
            self.vis.update_geometry(self.cloud)
            self.vis.poll_events()
            self.vis.update_renderer()
        return False

    def _prev(self,
        vis:o3d.visualization.Visualizer
        ) -> bool:
        """
        Callback to go to the previous iteration in the 3D viewer.
    
        Args:
            vis (open3d.visualization.Visualizer): The Open3D visualizer instance triggering the callback.
    
        Returns:
            bool: False to indicate that the visualizer should continue running.
        """
        if self.iteration > 0:
            self.iteration -= 1
            self._update_geometry()
            self.vis.update_geometry(self.cloud)
            self.vis.poll_events()
            self.vis.update_renderer()
        return False

    def _random_color_ID(self,
        ID:int
        ) -> np.ndarray:
        """
        Generates a reproducible random RGB color based on the given ID.
    
        Args:
            ID (int): A unique identifier used to seed the random generator.
    
        Returns:
            np.ndarray: An array of 3 floats representing an RGB color in [0, 1].
        """
        np.random.seed(ID)
        return np.random.rand(3)

    def _print_progress_bar(self,
        iteration:int,
        total:int,
        length:int = 40
        ) -> None:
        """
        Display the progress bar of the simulation on the terminal.
        
        Args:
            iteration (int): Iteration to show.
            total (int): Number of iterations.
            length (int, optional): Lentgh of the progress bar displayed in the terminal.
        """
        percent = iteration / (total - 1) if total > 1 else 1
        filled = int(length * percent)
        bar = '█' * filled + '-' * (length - filled)
        text = f"Progress: |{bar}| {iteration + 1}/{total}"
    
        sys.stdout.write('\r' + text.ljust(80))
        sys.stdout.flush()

# =============================================================================
# CLASS FragmentDetector
# =============================================================================
class FragmentDetector:
    def __init__(self,
        source:Union[str, List[str]],
        trash_treshold:int = 10):
        """
        This class loads AGDD files from a directory or list of file paths,
        and provides methods to detect, visualize, and analyze fragments.
        
        Args:
            source (str | list[str]): Path to a directory of AGDD files or a list of file paths.
            trash_treshold (int, optional): Minimum number of points to consider a fragment.
    
        Raises:
            ValueError: If `source` is not a string or a list of paths.
        """
        self._trash_treshold = trash_treshold

        if isinstance(source, str):
            self._directory = source
            self._agdd_files = self._load_agdd()
        elif isinstance(source, list):
            self._directory = os.path.commonpath(source)
            self._agdd_files = sorted(source)
        else:
            raise ValueError("Source must be a directory path or a list of AGDD files.")
        
        self._fragments = []
        self._iterations_nb = len(self._agdd_files)
        self._seed = random.randint(0, 2**32 - 1)
        self._2d_plot_lim = []
        
# === Core methods ===
    def build_fragments(self) -> None:
        """
        Detects fragments from the loaded AGDD files and stores them internally.
        """
        for iteration, agdd_file in enumerate(self._agdd_files):
            new_fragments = Fragment.detect_from_agdd(self, agdd_file, iteration, self._trash_treshold)
            self._fragments.append(new_fragments)
    
    def get_fragments_by_ID(self, 
        ID:int
        ) -> List['Fragment']:
        """
        Get all fragments that have the specified ID across all iterations.

        Args:
            ID (int): Fragment ID.

        Returns:
            List[Fragment]: List of fragments with the given ID.
        """
        list_ = []
        for fragment_list in self._fragments:
            for fragment in fragment_list:
                if fragment.ID == ID:
                    list_.append(fragment)
        return list_
    
    def get_fragments_by_iteration(self,
        iteration:int
        ) -> List['Fragment']:
        """
        Return fragments for a specific iteration.

        Args:
            iteration (int): Iteration index (must be >= 0).

        Returns:
            List[Fragment]: List of fragments for the given iteration.

        Raises:
            ValueError: If iteration is negative.
            IndexError: If iteration is out of range.
        """
        if iteration < 0:
            raise ValueError("Iteration must be non-negative.")
        if iteration >= len(self._fragments):
            raise IndexError(f"Iteration {iteration} is out of range (max is {len(self._fragments) - 1}).")
        return self._fragments[iteration]
    
    def save(self) -> None:
        """
        Save all fragment data to text files in a 'txt' subfolder.

        Raises:
            ValueError: If no fragments have been built.
        """
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

# === Visualization ===
    def plot2D(self,
        iteration:int,
        show_trash:bool = False,
        show_bonds:bool = False,
        auto_close:bool = False,
        dot_size:int = 4,
        rigid_lim:bool = True,
        border_coef:float = 1.5,
        save:bool = True,
        save_format:str = "png"
        ) -> None:
        """
        Plot a 2D scatter plot of fragment nodes for a given iteration.

        Args:
            iteration (int): Iteration index to plot.
            show_trash (bool): Whether to show trash fragments (ID = -1).
            show_bonds (bool): Whether to show bonds between nodes.
            auto_close (bool): Whether to automatically close the plot window.
            dot_size (int): Size of the dots for each node.
            rigid_lim (bool): Whether to keep the axis limits fixed after first plot.
            border_coef (float): Expansion factor for plot borders.
            save (bool): Whether to save the plot as an image.
            save_format (str): File format to save image (e.g., "png", "pdf").

        Raises:
            ValueError: If no fragments were built before plotting.
        """
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
        xgrid:bool = True,
        interline:bool = True,
        auto_close:bool = False,
        save:bool = True,
        save_format:str = "png"
        ) -> None:
        """
        Plot a stack graph of fragment repartition over all iterations.

        Args:
            xgrid (bool): Wether to show vertical line marking iterations.
            interline (bool): Whether to show line between fragments repartition.
            auto_close (bool): Whether to automatically close the plot window.
            save (bool): Whether to save the plot as an image.
            save_format (str): File format to save image (e.g., "png", "pdf").

        Raises:
            ValueError: If no fragments were built before plotting.
        """
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
        auto_close:bool = False,
        save:bool = True,
        save_format:str = "png"
        ) -> None:
        """
        Plot a network graph of fragments.

        Args:
            auto_close (bool): Whether to automatically close the plot window.
            save (bool): Whether to save the plot as an image.
            save_format (str): File format to save image (e.g., "png", "pdf").

        Raises:
            ValueError: If no fragments were built before plotting.
        """
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
    
    def viewer3D(self) -> None:
        """
        Launch a 3D viewer to display detected fragments.
        
        Raises:
            ValueError: If no fragments have been built.
        """
        if not self.fragments:
            raise ValueError("No fragments built!")
            
        viewer = Viewer3D(self)
        viewer.run()

# === Internal helpers ===
    def _load_agdd(self) -> List[str]:
        """
        Load the AGDD files from the work directory.
    
        Returns:
            List[str]: List of AGDD paths.
        """
        agdd_files = glob.glob(os.path.join(self._directory, "*.agdd"))
        agdd_files.sort()
        return agdd_files
    
    def _get_fragments_by_iteration_safe(self,
        iteration:int
        ) -> List['Fragment']:
        """
        Return a list of fragments at a given iteration.
    
        Args:
            iteration (int): The iteration at which the fragments are detected.
    
        Returns:
            List['Fragment']: List of fragments detected at the given iteration.
        """
        if 0 <= iteration < len(self._fragments):
            return self._fragments[iteration]
        else:
            return []
    
    def _random_color_ID(self,
        ID:int
        ) -> tuple[float]:
        """
        Generates a reproducible random RGB color based on the given ID.
    
        Args:
            ID (int): A unique identifier used to seed the random generator.
    
        Returns:
            np.ndarray: An array of 3 floats representing an RGB color in [0, 1].
        """
        np.random.seed(self._seed + ID)
        r, g, b = np.random.rand(3)
        return (r, g, b)

# === Properties ===
    @property
    def directory(self) -> str:
        """
        str: Work directory from which AGDD files are loaded.
        """
        return self._directory
    
    @property
    def trash_treshold(self) -> int:
        """
        int: Minimum number of points for a valid fragment.
        """
        return self._trash_treshold
    
    @property
    def agdd_files(self) -> List[str]:
        """
        List[str]: List of AGDD file paths.
        """
        return self._agdd_files
    
    @property
    def fragments(self) -> List[List['Fragment']]:
        """
        List[List[Fragment]]: All detected fragments by iteration.
        """
        return self._fragments
    
    @property
    def iterations_nb(self) -> int:
        """
        int: Number of AGDD files processed.
        """
        return self._iterations_nb
    
# =============================================================================
# CLASS Fragment
# =============================================================================
class Fragment:
    def __init__(self,
        ID:int,
        nodeIDs:list,
        nodes:list,
        bondIDs:list,
        bonds:list,
        iteration:int,
        parentIDs:list,
        area:float):
        """
        Represents a fragment detected in the simulation.
        
        Args:
            ID (int): Unique identifier for the fragment.
            nodeIDs (np.ndarray[int]): Indices of the nodes that belong to this fragment.
            nodes (np.ndarray[float]): Node data (e.g., coordinates and radii).
            bondIDs (np.ndarray[int]): Indices of the bonds associated with this fragment.
            bonds (np.ndarray[int]): Bond data (e.g., node pairs forming bonds).
            iteration (int): Iteration number when this fragment was detected.
            parentIDs (List[int]): List of fragment IDs from the previous iteration that contributed to this one.
            area (float): Percentage of the total area covered by this fragment.
        """
        self._ID = ID
        self._nodeIDs = nodeIDs
        self._nodes = nodes
        self._bondIDs = bondIDs
        self._bonds = bonds
        self._iteration = iteration
        self._parentIDs = parentIDs
        self._area = area
        
# === Core methods ===
    @classmethod
    def detect_from_agdd(cls,
         detector:'FragmentDetector',
         agdd_file:str,
         iteration:int,
         trash_treshold:int
         )-> List['Fragment']:
        """
        Detects fragments from a given AGDD file.
    
        Args:
            detector (FragmentDetector): The main detector instance containing fragment and configuration data.
            agdd_file (str): Path to the AGDD file to process.
            iteration (int): Index of the current iteration (frame).
            trash_treshold (int): Minimum number of points required to consider a group as a valid fragment.
    
        Returns:
            List[Fragment]: A list of detected `Fragment` instances for the given iteration.
        """
        nodes_nb, bonds_nb, nodes_all, bonds_all = cls._read_agdd(agdd_file)
        connectivity = cls._compute_connectivity(nodes_nb, bonds_all)
        fragments = cls._identify_fragments(detector, connectivity, nodes_all, bonds_all, iteration, trash_treshold)
        return fragments
    
# === Internal helpers ===
    @staticmethod
    def _read_agdd(
        agdd_file:str
        ) -> tuple:
        """
        Extract data from an AGDD file.
    
        Args:
            agdd_file (str): Path to the AGDD file.
    
        Returns:
            tuple: A tuple containing:
                - nodes_nb (int): Number of nodes.
                - bonds_nb (int): Number of bonds.
                - nodes_all (np.ndarray): Array of node coordinates.
                - bonds_all (np.ndarray): Array of bond connections.
        """
        file = open(agdd_file, "r")
        lines = file.readlines()
        file.close()
        
        nodes_nb = int(lines[0])
        bonds_nb = int(lines[nodes_nb + 1])
        
        nodes_all = np.array([list(map(float, line.split())) for line in lines[1:nodes_nb + 1]])
        bonds_all = np.array([list(map(int, line.split())) for line in lines[nodes_nb + 2:nodes_nb + 2 + bonds_nb]])
        return nodes_nb, bonds_nb, nodes_all, bonds_all
    
    @staticmethod
    def _compute_connectivity(
        nodes_nb:int,
        bonds_all:list
        ) -> csr_matrix:
        """
        Creates a symmetric node connectivity matrix from bond data.
    
        Args:
            nodes_nb (int): Total number of nodes.
            bonds_all (List[List[int]]): List of bonds, each represented as a pair of connected node indices.
    
        Returns:
            csr_matrix: A sparse adjacency matrix (Compressed Sparse Row format) representing node connectivity.
        """
        connectivity = lil_matrix((nodes_nb, nodes_nb), dtype=int)
        for i, j in bonds_all:
            connectivity[i, j] = 1
            connectivity[j, i] = 1
        return connectivity.tocsr()
    
    def _identify_fragments(
        detector:'FragmentDetector',
        connectivity:csr_matrix,
        nodes_all:list,
        bonds_all:list,
        iteration:int,
        trash_treshold:int
        ) -> List['Fragment']:
        """
        Identifies fragments from a connectivity graph of nodes and bonds.
    
        Args:
            detector (FragmentDetector): The main detector instance containing prior fragment data.
            connectivity (csr_matrix): Sparse adjacency matrix representing node connectivity.
            nodes_all (List[List[float]]): All nodes, including coordinates and radii.
            bonds_all (List[List[int]]): All bonds, defined by pairs of connected node indices.
            iteration (int): The current iteration index.
            trash_treshold (int): Minimum number of nodes required for a fragment to be considered valid.
    
        Returns:
            List[Fragment]: A list of `Fragment` objects, including valid fragments and possibly a trash fragment.
        """
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

# === Properties ===
    @property
    def ID(self) -> int:
        """
        int: Unique identifier of the fragment.
        """
        return self._ID
    
    @property
    def nodeIDs(self) -> np.ndarray[int]:
        """
        np.ndarray[int]: Indices of the nodes that belong to the fragment.
        """
        return self._nodeIDs
    
    @property
    def nodes(self) -> np.ndarray[float]:
        """
        np.ndarray[float]: Array containing the coordinates and radii of the fragment's nodes.
        """
        return self._nodes
    
    @property
    def bondIDs(self) -> np.ndarray[int]:
        """
        np.ndarray[int]: Indices of the bonds associated with the fragment.
        """
        return self._bondIDs
    
    @property
    def bonds(self) -> np.ndarray[int]:
        """
        np.ndarray[int]: Array of node connections (bonds) within the fragment.
        """
        return self._bonds
    
    @property
    def iteration(self) -> int:
        """
        int: Iteration number at which the fragment is detected.
        """
        return self._iteration
    
    @property
    def parentIDs(self) -> List[int]:
        """
        List[int]: IDs of parent fragments from the previous iteration.
        """
        return self._parentIDs
    
    @property
    def area(self) -> float:
        """
        float: Percentage of the total area occupied by the fragment.
        """
        return self._area
    
    

    
    
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





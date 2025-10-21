# =============================================================================
# LIB
# =============================================================================
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from collections import defaultdict
import matplotlib.pyplot as plt
import pyvista as pv
import random
from pathlib import Path
import os

from .classes import Node, Bond, Fragment

# =============================================================================
# GLOBAL VARIABLE
# =============================================================================
SEED = random.randint(0, 2**32 - 1)

# =============================================================================
# GENERAL FUNCTIONS
# =============================================================================
def load_files(
        directory: str,
        extension: str = "txt"
        ) -> list[str]:
    """
    Load files from a specified directory.

    Parameters
    ----------
    directory : str
        Search directory.
    extension : str, optional
        Extension of the files to list. The default is "txt".

    Returns
    -------
    list[str]
        List of files.
    """
    path = Path(directory)
    
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"The directory '{directory}' does not exist.")
    
    files = sorted([str(p) for p in path.glob(f"*{extension}")])
    
    if not files:
        raise ValueError(f"No files with extension '{extension}' found in '{directory}'.")
    
    return files
    
def detect_fragment(
        files: list[str],
        trash_threshold: int = 10
        ) -> dict[int, list[Fragment]]:
    """
    Compute fragments from a list of files.

    Parameters
    ----------
    files : list[str]
        List of files to analyze.
    trash_threshold : int, optional
        Minimum fragment size, all fragment of size below the threshold will be considered trash. The default is 10.

    Returns
    -------
    dict[int, list[Fragment]]
        Dictionnary of fragments.
    """
    fragdict = {}
    prevfrag = []
    next_id = 0
    
    for it, file in enumerate(files):
        # read file
        n_all, b_all = _read_file(file)
        
        # connectivity matrix
        connectivity = _compute_connectivity(len(n_all), b_all)
        
        # compute fragments
        num_frag, n_label_all = connected_components(connectivity, directed=False)
        
        # reassign small fragments to trash (-1)
        label_counts = np.bincount(n_label_all)
        n_label_all = np.array([n_label if label_counts[n_label] >= trash_threshold else -1 
                                for n_label in n_label_all])
        
        # group nodes and bonds per fragment
        frag_n = defaultdict(list)
        frag_b = defaultdict(list)
        for n_i, n_label in enumerate(n_label_all):
            frag_n[n_label].append(n_all[n_i])
        for b in b_all.values():
            n1_label = n_label_all[b.n1.i]
            n2_label = n_label_all[b.n2.i]
            if n1_label == n2_label:
                frag_b[n1_label].append(b)
        
        # build fragment objects
        fraglist = []
        frag_i = set()
        
        for label in sorted(set(n_label_all)):
            n_list = frag_n.get(label, [])
            b_list = frag_b.get(label, [])
            ratio = 100. * len(n_list) / len(n_all)
            
            # compute ancestors
            ancestor = []
            n_i = set(n.i for n in n_list)
            for prev_frag in prevfrag:
                prev_n_i = set(n.i for n in prev_frag._nodes)
                if len(n_i & prev_n_i) > 0:
                    ancestor.append(prev_frag._i)
            
            # determine fragment ID
            if it == 0:
                frag_id = -1 if label == -1 else next_id
                if label != -1:
                    next_id += 1
            else:
                if label == -1:
                    frag_id = -1
                elif ancestor:
                    main_ancestor = ancestor[0]
                    if main_ancestor in frag_i:
                        frag_id = next_id
                        next_id += 1
                    else:
                        frag_id = main_ancestor
                        frag_i.add(main_ancestor)
                else:
                    frag_id = next_id
                    next_id += 1
            
            fraglist.append(Fragment(frag_id, n_list, b_list, it, ratio, ancestor))
        
        fragdict[it] = fraglist
        prevfrag = fraglist
    
    return fragdict

def get_fragments_by_id(
        fragments: dict[int, list[Fragment]] | list[Fragment]
        ) -> dict[int, list[Fragment]]:
    """
    Return a dictionnary of fragments sorted by id.

    Parameters
    ----------
    fragments : dict[int, list[Fragment]] | list[Fragment]
        Fragment objects.

    Returns
    -------
    dict[int, list[Fragment]]
        Dictionnary of fragments sorted by id.
    """
    fragdict = {}
    
    if isinstance(fragments, dict):
        fraglist = [f for flist in fragments.values() for f in flist]
    elif isinstance(fragments, list):
        fraglist = fragments.copy()
    else:
        raise TypeError("The fragments must be a dictionnary or a list of Fragment objects")
    
    for f in fraglist:
        try:
            fragdict[f.i].append(f)
        except KeyError:
            fragdict[f.i] = [f]
                    
    return fragdict

def get_fragments_by_it(
        fragments: dict[int, list[Fragment]] | list[Fragment]
        ) -> dict[int, list[Fragment]]:
    """
    Return a dictionnary of fragments sorted by iteration.

    Parameters
    ----------
    fragments : dict[int, list[Fragment]] | list[Fragment]
        Fragment objects.

    Returns
    -------
    dict[int, list[Fragment]]
        Dictionnary of fragments sorted by iteration.
    """
    fragdict = {}
    
    if isinstance(fragments, dict):
        fraglist = [f for flist in fragments.values() for f in flist]
    elif isinstance(fragments, list):
        fraglist = fragments.copy()
    else:
        raise TypeError("The fragments must be a dictionnary or a list of Fragment objects")
    
    for f in fraglist:
        try:
            fragdict[f.it].append(f)
        except KeyError:
            fragdict[f.it] = [f]
                    
    return fragdict

# =============================================================================
# INTERNAL HELPERS
# =============================================================================
def _read_file(
        file: str
        ) -> tuple[dict[int, Node], dict[int, Bond]]:
    """
    Internal helper to read a file.

    Parameters
    ----------
    file : str
        File to read (containing nodes positions and bond connectivity).

    Returns
    -------
    (tuple[dict[int, Node], dict[int, Bond]])
        Dictionaries of nodes and bonds.
    """
    with open(file, "r") as f:
        lines = f.readlines()
    
    num_n = int(lines[0])
    num_b = int(lines[num_n + 1])
    
    n_all = {}
    for i, line in enumerate(lines[1:num_n + 1]):
        x, y, z, *_ = list(map(float, line.split()))
        n_all[i] = Node(i, x, y, z)

    b_all = {}
    for i, line in enumerate(lines[num_n + 2:num_n + 2 + num_b]):
        n1, n2 = list(map(int, line.split()))
        b_all[i] = Bond(i, n_all[n1], n_all[n2])
    
    return n_all, b_all

def _compute_connectivity(
        num_n: int,
        bonds_all: dict[int, Bond]
        ) -> csr_matrix:
    """
    Compute a light connectivity matrix.

    Parameters
    ----------
    num_n : int
        Number of nodes.
    bonds_all : dict[int, Bond]
        Dictionnary of bonds that connect the nodes.

    Returns
    -------
    csr_matrix
        Matrix of connectivity.
    """
    connectivity = lil_matrix((num_n, num_n), dtype=int)
    for bond in bonds_all.values():
        connectivity[bond.n1.i, bond.n2.i] = 1
        connectivity[bond.n2.i, bond.n1.i] = 1
    return connectivity.tocsr()

def _random_color_id(
        i: int,
        seed: int = SEED
        ) -> tuple[float]:
    """
    Generate a random color within a colormap randomly initialized.

    Parameters
    ----------
    i : int
        Indice of the element to color.
    seed : int, optional
        Seed of the colormap. The default is SEED.

    Returns
    -------
    tuple[float]
        R, G, B color values.
    """
    np.random.seed(seed + i)
    r, g, b = np.random.rand(3)
    return (r, g, b)

# =============================================================================
# VISUALIZATION TOOLKIT
# =============================================================================
def pvplot(
        fraglist: list[Fragment], 
        xlim: list[int] | None = None,
        ylim: list[int] | None = None,
        zlim: list[int] | None = None,
        show_node: bool = True,
        show_bond: bool = True,
        show_axe: bool = False,
        auto_close: bool = False,
        camera_position: list[tuple] | None = None,
        save: bool = True,
        filename: str = "plot",
        save_dir: str = os.path.join(os.getcwd(), "img")
        ) -> None:
    """
    Plot fragments.

    Parameters
    ----------
    fraglist : list[Fragment]
        List of Fragment objects.
    xlim : list[int] | None, optional
        Set plot x-limits. The default is None.
    ylim : list[int] | None, optional
        Set plot y-limits. The default is None.
    zlim : list[int] | None, optional
        Set plot z-limits. The default is None.
    show_node : bool, optional
        Show nodes. The default is True.
    show_bond : bool, optional
        Show bonds. The default is True.
    show_axe : bool, optional
        Show axes. The default is False.
    auto_close : bool, optional
        Automatically close window plot. The default is False.
    camera_position : list[tuple] | None, optional
        Set the camera position (see https://docs.pyvista.org/api/plotting/_autosummary/pyvista.cameraposition#pyvista.CameraPosition). The default is 'xy'.
    save : bool, optional
        Save option. The default is True.
    filename : str
        Filename of saved image.
    save_dir : str, optional
        Save directory. The default is "/img" in current directory.
    """
    plotter = pv.Plotter(off_screen=auto_close)
    
    for frag in fraglist:
        # fragment color
        color = _random_color_id(frag.i)
        
        # build points
        points = []
        id_to_local = {}
        for i, n in enumerate(frag.nodes):
            # node coordinates
            points.append((n.x, n.y, n.z))
            
            # mapping: global node id -> local index (0..N-1)
            id_to_local[n.i] = i
            
        points = np.array(points, dtype=float)
        
        # plot nodes
        n_pdata = pv.PolyData(points)
        if show_node: plotter.add_mesh(n_pdata, render_points_as_spheres=True, point_size=8, color=color)
        
        # plot fragment indice
        plotter.add_point_labels([frag.centroid.x, frag.centroid.y, frag.centroid.z], 
                                 [frag.i], font_size=20, text_color=color, 
                                 always_visible=True, shape_opacity=1.0, 
                                 shape_color="white")
        
        # build bonds
        bonds = []
        for b in frag.bonds:
            n1 = b.n1.i
            n2 = b.n2.i
            bonds.append((id_to_local[n1], id_to_local[n2]))
        
        # build lines
        if bonds:
            lines = np.hstack([[2, i, j] for (i, j) in bonds])
            
            b_pdata = pv.PolyData()
            b_pdata.points = points
            b_pdata.lines = lines

            # plot bonds
            if show_bond: plotter.add_mesh(b_pdata, color=color, line_width=2)

    if show_axe: 
        def auto_tick_format(r): return "%.0f" if r > 100 else ("%.2f" if r > 1 else "%.4f")
        
        xmin, ymin, zmin = points.min(axis=0)
        xmax, ymax, zmax = points.max(axis=0)
    
        plotter.show_bounds(
            grid = 'front', 
            font_size = 10,
            fmt=auto_tick_format(np.max([xmax - xmin, ymax - ymin, zmax - zmin])))
    
    if camera_position == None:
        plotter.camera_position = 'xy'
    else:
        plotter.camera_position = camera_position
        
    plotter.show()
    
    if save:
        filepath = os.path.join(save_dir, filename)
        plotter.screenshot(filepath)
        print(f"Plot saved to {filepath}")
    
def stackplot(
        fraglist: list[Fragment], 
        xgrid: bool = True,
        contour: bool = True,
        auto_close: bool = False,
        save: bool = True,
        filename: str = "stackplot", 
        save_dir: str = os.path.join(os.getcwd(), "img")
        ) -> None:
    """
    Stackplot of the evolution of the fragment repartition.

    Parameters
    ----------
    fragments : list[Fragment]
        Fragment list.
    xgrid : bool, optional
        Show iteration separator. The default is True.
    contour : bool, optional
        Contour line of ech fragment. The default is True.
    auto_close : bool, optional
        Automatically close window plot. The default is False.
    save : bool, optional
        Save option. The default is True.
    filename : str
        Filename of saved image.
    save_dir : str, optional
        Save directory. The default is "/img" in current directory.
    """
    if len(fraglist) == 0:
        raise ValueError("The fragment list is empty")
    
    if not all(isinstance(f, Fragment) for f in fraglist):
        raise TypeError("The fragment list must contain Fragment objects")
    
    # get iterations and indices
    ids = sorted(get_fragments_by_id(fraglist).keys())
    its = sorted(get_fragments_by_it(fraglist).keys())

    # data matrix
    data = np.zeros((len(ids), len(its)))
    id_to_row = {fid: i for i, fid in enumerate(ids)}
    it_to_col = {fit: it for it, fit in enumerate(its)}

    for f in fraglist:
        if f.i in id_to_row and f.it in it_to_col:
            row = id_to_row[f.i]
            col = it_to_col[f.it]
            data[row, col] = f.ratio
    cumulative_data = np.cumsum(data, axis=0)
    
    # colors
    colors = [_random_color_id(fid) for fid in ids]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.stackplot(its, data, labels=ids, colors=colors, alpha=0.8)
    if xgrid:
        for it in its: ax.axvline(x=it, color='k', linestyle='--', linewidth=0.5)
    if contour:
        for line in cumulative_data[:-1]: ax.plot(its, line, color="k", linewidth=1)
    
    ax.legend(loc='center left', framealpha=1., bbox_to_anchor=(1., 0.5))
    ax.set_xlim(0, np.max(its))
    ax.set_ylim(0, 100)
    ax.set_xlabel("Iteration", fontsize=10)
    ax.set_ylabel("Repartition \ %", fontsize=10)
    plt.tight_layout()
    
    if save:
        os.makedirs(save_dir, exist_ok=True)
        filepath = os.path.join(save_dir, filename)
        plt.savefig(filepath)
        print(f"Plot saved to {filepath}")
    
    if auto_close:
        plt.close(fig)
    else:
        plt.show()

def histplot(
        fraglist: list[Fragment], 
        auto_close: bool = False,
        save: bool = True,
        filename: str = "histplot", 
        save_dir: str = os.path.join(os.getcwd(), "img")
        ) -> None:
    """
    Plot of the distribution of the fragment properties.

    Parameters
    ----------
    fragments : list[Fragment]
        Fragments list.
    auto_close : bool, optional
        Automatically close window plot. The default is False.
    save : bool, optional
        Save option. The default is True.
    filename : str
        Filename of saved image.
    save_dir : str, optional
        Save directory. The default is "/img" in current directory.
    """
    if len(fraglist) == 0:
        raise ValueError("The fragment list is empty")
    
    if not all(isinstance(f, Fragment) for f in fraglist):
        raise TypeError("The fragment list must contain Fragment objects")
        
    if not all(f.dim in (2, 3) for f in fraglist):
        raise ValueError("The fragments must have the same dimension")
    
    # get indices
    ids = sorted(get_fragments_by_id(fraglist).keys())

    # data matrix
    data = np.zeros((len(ids), 6))
    id_to_row = {fid: i for i, fid in enumerate(ids)}

    for f in fraglist:
        if f.i in id_to_row:
            row = id_to_row[f.i]
            data[row, 0] = len(f.nodes)
            data[row, 1] = len(f.bonds)
            data[row, 2] = f.ratio
            data[row, 3] = f.sphericity
            data[row, 4] = f.area
            data[row, 5] = f.volume if f.dim == 3 else f.perimeter
    
    fig, axs = plt.subplots(2, 3, figsize=(10, 8))
    
    # nodes
    counts, bins, patches = axs[0, 0].hist(data[:, 0], bins='auto')
    axs[0, 0].clear()
    axs[0, 0].bar(bins[:-1], counts / counts.sum() * 100, width=np.diff(bins), color='C0', edgecolor='black', alpha=0.7)
    axs[0, 0].set_xlabel("Num. of nodes")
    axs[0, 0].set_ylabel("Frag. percentage \ %")
    
    # bonds
    counts, bins, patches = axs[0, 1].hist(data[:, 1], bins='auto')
    axs[0, 1].clear()
    axs[0, 1].bar(bins[:-1], counts / counts.sum() * 100, width=np.diff(bins), color='C1', edgecolor='black', alpha=0.7)
    axs[0, 1].set_xlabel("Num. of bonds")
    axs[0, 1].set_ylabel("Frag. percentage \ %")
    
    # ratio
    counts, bins, patches = axs[0, 2].hist(data[:, 2], bins='auto')
    axs[0, 2].clear()
    axs[0, 2].bar(bins[:-1], counts / counts.sum() * 100, width=np.diff(bins), color='C2', edgecolor='black', alpha=0.7)
    axs[0, 2].set_xlabel("Space ratio \ %")
    axs[0, 2].set_ylabel("Frag. percentage \ %")
    
    # sphericity
    counts, bins, patches = axs[1, 0].hist(data[:, 3], bins='auto')
    axs[1, 0].clear()
    axs[1, 0].bar(bins[:-1], counts / counts.sum() * 100, width=np.diff(bins), color='C3', edgecolor='black', alpha=0.7)
    axs[1, 0].set_xlabel("Sphericity")
    axs[1, 0].set_ylabel("Frag. percentage \ %")
    
    # area
    counts, bins, patches = axs[1, 1].hist(data[:, 4], bins='auto')
    axs[1, 1].clear()
    axs[1, 1].bar(bins[:-1], counts / counts.sum() * 100, width=np.diff(bins), color='C4', edgecolor='black', alpha=0.7)
    axs[1, 1].set_xlabel("Area \ $m^2$")
    axs[1, 1].set_ylabel("Frag. percentage \ %")
    
    # volume / perimeter
    counts, bins, patches = axs[1, 2].hist(data[:, 5], bins='auto')
    axs[1, 2].clear()
    axs[1, 2].bar(bins[:-1], counts / counts.sum() * 100, width=np.diff(bins), color='C5', edgecolor='black', alpha=0.7)
    axs[1, 2].set_xlabel("Volume \ $m^3$") if f.dim == 3 else axs[1, 2].set_xlabel("Perimeter \ m")
    axs[1, 2].set_ylabel("Frag. percentage \ %")
    
    plt.tight_layout()
    plt.show()

    if save:
        os.makedirs(save_dir, exist_ok=True)
        filepath = os.path.join(save_dir, filename)
        plt.savefig(filepath)
        print(f"Plot saved to {filepath}")
    
    if auto_close:
        plt.close(fig)
    else:
        plt.show()

def stats(
        fraglist: list[Fragment], 
        save: bool = True,
        filename: str = "stats", 
        save_dir: str = os.path.join(os.getcwd())
        ) -> None:
    """
    Retun the statistics of the fragment list.
    
    Parameters
    ----------
    fragments : list[Fragment]
        Fragments list.
    save : bool, optional
        Save option. The default is True.
    filename : str
        Filename of saved stats.
    save_dir : str, optional
        Save directory. The default is the current directory.
    """
    if len(fraglist) == 0:
        raise ValueError("The fragment list is empty")
    
    if not all(isinstance(f, Fragment) for f in fraglist):
        raise TypeError("The fragment list must contain Fragment objects")
    
    if not all(f.dim in (2, 3) for f in fraglist):
        raise ValueError("The fragments must have the same dimension")
    
    # get indices
    ids = sorted(get_fragments_by_id(fraglist).keys())
    
    # data matrix
    data = np.zeros((len(ids), 6))
    id_to_row = {fid: i for i, fid in enumerate(ids)}
    
    for f in fraglist:
        if f.i in id_to_row:
            row = id_to_row[f.i]
            data[row, 0] = len(f.nodes)
            data[row, 1] = len(f.bonds)
            data[row, 2] = f.ratio
            data[row, 3] = f.sphericity
            data[row, 4] = f.area
            data[row, 5] = f.volume if f.dim == 3 else f.perimeter
    
    # property names
    prop_names = ["Num. of nodes", "Num. of bonds", "Space ratio", 
                  "Sphericity", "Area"]
    prop_names += ["Perimeter" if f.dim == 2 else "Volume"]
    
    # compute statistics
    stats_dict = {}
    for i, name in enumerate(prop_names):
        vals = data[:, i]
        stats_dict[name] = {
            "mean": np.mean(vals),
            "median": np.median(vals),
            "min": np.min(vals),
            "max": np.max(vals),
            "std": np.std(vals)
        }
    
    # print stats
    for name, s in stats_dict.items():
        print(f"--- {name} ---")
        for k, v in s.items():
            print(f"{k.capitalize():>6}: {v:.4f}")
        print()
    
    # save to file if requested
    if save:
        os.makedirs(save_dir, exist_ok=True)
        if not filename.lower().endswith(".txt"):
            filename += ".txt"
        filepath = os.path.join(save_dir, filename)
        with open(filepath, "w") as f:
            for name, s in stats_dict.items():
                f.write(f"--- {name} ---\n")
                for k, v in s.items():
                    f.write(f"{k.capitalize():>6}: {v:.4f}\n")
                f.write("\n")
        print(f"Statistics saved to {filepath}")
            
# =============================================================================
# DEBUGGING
# =============================================================================
if __name__ == "__main__":
    pass








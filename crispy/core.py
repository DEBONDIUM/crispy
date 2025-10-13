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
def load_files(directory: str, extension: str = "txt") -> list[str]:
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
    files = [str(p) for p in Path(directory).glob(f"*{extension}")]
    files.sort()
    return files
    
def detect_fragment(files: list[str], trash_threshold: int = 10) -> dict[int, list[Fragment]]:
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
    fragments_history = {}
    prev_fragments = []
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
            for prev_frag in prev_fragments:
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
        
        fragments_history[it] = fraglist
        prev_fragments = fraglist
    
    return fragments_history

def get_fragments_by_id(fragments: dict[int, list[Fragment]] | list[Fragment]) -> dict[int, list[Fragment]]:
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

def get_fragments_by_it(fragments: dict[int, list[Fragment]] | list[Fragment]) -> dict[int, list[Fragment]]:
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
def _read_file(file: str) -> tuple[dict[int, Node], dict[int, Bond]]:
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

def _compute_connectivity(num_n: int, bonds_all: dict[int, Bond]) -> csr_matrix:
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


def _random_color_id(i: int, seed: int = SEED) -> tuple[float]:
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
def plot(fraglist: list[Fragment], 
           title: str = "",
           xlim: list[int] | None = None,
           ylim: list[int] | None = None,
           zlim: list[int] | None = None,
           show_node: bool = True,
           show_bond: bool = True,
           show_axe: bool = False,
           auto_close: bool = False,
           save: bool = True,
           filename: str = "plot",
           save_dir: str = os.path.join(os.getcwd(), "img")) -> None:
    """
    Plot fragments.

    Parameters
    ----------
    fraglist : list[Fragment]
        List of Fragment objects.
    title : str, optional
        Plot title. The default is "".
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
    save : bool, optional
        Save option. The default is True.
    filename : str
        Filename of saved image.
    save_dir : str, optional
        Save directory. The default is os.path.join(os.getcwd(), "img").
    """
    plotter = pv.Plotter(title=title, off_screen=auto_close)
    
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
    
    plotter.camera_position = [
        (0, 0, 10),   # camera location (x, y, z)
        (0, 0, 0),    # focal point (center of view)
        (0, 1, 0)     # up vector (keeps Y upward)
    ]
    
    plotter.show(cpos='xy')
    
    if save:
        os.makedirs(save_dir, exist_ok = True)
        plotter.screenshot(os.path.join(save_dir, filename))
    
def stackplot(fragments: dict[int, list[Fragment]] | list[Fragment], 
           title: str = "",
           xgrid: bool = True,
           contour: bool = True,
           auto_close: bool = False,
           save: bool = True,
           filename: str = "stackplot", 
           save_dir: str = os.path.join(os.getcwd(), "img")) -> None:
    """
    Stackplot of the evolution of the fragment repartition.

    Parameters
    ----------
    fragments : dict[int, list[Fragment]] | list[Fragment]
        Fragment objetcs.
    title : str, optional
        Figure title. The default is "".
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
        Save directory. The default is os.path.join(os.getcwd(), "img").
    """
    if isinstance(fragments, dict):
        fraglist = [f for flist in fragments.values() for f in flist]
    elif isinstance(fragments, list):
        fraglist = fragments.copy()
    else:
        raise TypeError("The fragments must be a dictionnary or a list of Fragment objects")
    
    # get iterations and indices
    ids = sorted(get_fragments_by_id(fraglist).keys())
    its = sorted(get_fragments_by_it(fraglist).keys())

    # data matrix
    data = np.zeros((len(ids), len(its)))
    id_to_row = {fid: i for i, fid in enumerate(ids)}

    for f in fraglist:
        if f.i in id_to_row:
            row = id_to_row[f.i]
            data[row, f.it] = f.ratio
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
    fig.title = title
    plt.tight_layout()
    
    if save:
        os.makedirs(save_dir, exist_ok = True)
        plt.savefig(os.path.join(save_dir, filename))
    
    if auto_close:
        plt.close(fig)
    else:
        plt.show()
            
# =============================================================================
# DEBUGGING
# =============================================================================
if __name__ == "__main__":
    pass








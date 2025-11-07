import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import networkx as nx
from particle import Particle
import hipopy.hipopy as hp
import argparse
import awkward as ak
import numpy as np
import vector

import vector
import numpy as np

# Set font sizes
plt.rc('font', size=25) #controls default text size
plt.rc('axes', titlesize=50) #fontsize of the title
plt.rc('axes', labelsize=50) #fontsize of the x and y labels
plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
plt.rc('legend', fontsize=25) #fontsize of the legend

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True

def boost_event_dict_to_gammaN(event_dict, beam_lepton, scattered_lepton, target_mass=0.938):
    """
    Boost all particles in event_dict to the gamma*-N CM frame using vector package.
    
    Inputs:
        event_dict: dict with keys 'E','px','py','pz'
        beam_lepton: 4-vector [E, px, py, pz] of incoming electron
        scattered_lepton: 4-vector [E, px, py, pz] of scattered electron
        target_mass: proton mass in GeV
    Returns:
        boosted_event_dict: same keys but in gamma*N CM frame
        boost_vector: vector object for the boost
    """
    # --- Create 4-vectors ---
    particles = vector.arr({
        "E": event_dict["E"],
        "px": event_dict["px"],
        "py": event_dict["py"],
        "pz": event_dict["pz"]
    })

    k = vector.obj(E=beam_lepton[0], px=beam_lepton[1], py=beam_lepton[2], pz=beam_lepton[3])
    kprime = vector.obj(E=scattered_lepton[0], px=scattered_lepton[1], py=scattered_lepton[2], pz=scattered_lepton[3])
    q = k - kprime

    # Target nucleon
    pN = vector.obj(E=target_mass, px=0, py=0, pz=0)

    # Gamma*-N system 4-vector
    P_sys = q + pN

    # Boost all particles to CM frame
    boosted_particles = particles.boostCM_of_p4(P_sys)

    # Build boosted_event_dict
    boosted_event_dict = {
        "E": boosted_particles.E,
        "px": boosted_particles.px,
        "py": boosted_particles.py,
        "pz": boosted_particles.pz,
    }

    return boosted_event_dict, P_sys


def plot_decay_chain(event_dict,bank_prefix="MC::Lund_",event_idx=0,layout="hierarchy"):
    """
    Plots the decay chain of a particle physics event with PID-based node styles.

    Parameters
    ----------
    event_dict : dict
        Must contain:
            - bank_prefix+'pid'      : list or array of particle IDs
            - bank_prefix+'index'    : list or array of 1-based indices (1..nparticles)
            - bank_prefix+'parent'   : list or array of parent indices (1..nparticles)
            - bank_prefix+'daughter' : list or array of daughter indices (1..nparticles)
                          or 0 / None / -1 for root particle.
    """
    pids = event_dict[bank_prefix+'pid']
    indices = event_dict[bank_prefix+'index']
    parents = event_dict[bank_prefix+'parent']
    daughters = event_dict[bank_prefix+'daughter']

    G = nx.DiGraph()

    for idx, pid in zip(indices, pids):
        G.add_node(idx, pid=pid)

    for idx, parent_idx in zip(indices, parents):
        if parent_idx not in [None, 0, -1]:
            G.add_edge(parent_idx, idx)

    if layout != "mothers_only":
        for idx, daughter_idx in zip(indices, daughters):
            if daughter_idx not in [None, 0, -1]:
                G.add_edge(idx, daughter_idx)

    # --- Layout ---
    if layout == "mothers_only":
        # ---- Hierarchical layout ----
        def hierarchy_pos(G, root=None, width=0.6, vert_gap=0.3, vert_loc=0, xcenter=0.5):
            if not nx.is_tree(G):
                raise TypeError("Cannot use hierarchy_pos on a graph that is not a tree")

            if root is None:
                roots = [n for n, d in G.in_degree() if d == 0]
                if not roots:
                    raise ValueError("No root found (no node with in_degree=0)")
                root = roots[0]

            def _hierarchy_pos(G, root, left, right, vert_loc, vert_gap, pos):
                children = list(G.successors(root))
                pos[root] = ((left + right) / 2, vert_loc)
                if children:
                    dx = (right - left) / len(children)
                    nextx = left
                    for child in children:
                        nextx += dx
                        pos = _hierarchy_pos(G, child, nextx - dx, nextx,
                                            vert_loc - vert_gap, vert_gap, pos)
                return pos

            return _hierarchy_pos(G, root, 0, width, vert_loc, vert_gap, {})

        # try:
        #     pos = hierarchy_pos(G)
        # except Exception:
        #     pos = nx.spring_layout(G, seed=42)

        # For directed decay chain
        components = list(nx.weakly_connected_components(G))

        subgraph_positions = {}
        x_offset = 0.0  # horizontal offset for placing components

        for comp in components:
            H = G.subgraph(comp)
            pos = hierarchy_pos(H, width=2.5)
            # Get max x-coordinate in this component
            max_x = max(x for x, y in pos.values())
            # Shift component
            pos = {n: (x + x_offset, y) for n, (x, y) in pos.items()}
            subgraph_positions.update(pos)
            # Add spacing for next component
            x_offset += max_x + 1.0  # 1.0 unit gap

        pos = subgraph_positions

    if layout == "rapidities":

        # --- Regions ---
        q_dir = 1.0  # assume +z if not provided

        rapidities = []
        for i in range(len(event_dict["pid"])):
            E = event_dict["E"][i]
            pz = event_dict["pz"][i]
            try:
                y = 0.5 * np.log((E + pz*q_dir) / (E - pz*q_dir))
            except ZeroDivisionError:
                y = np.nan
            rapidities.append(y)

        # --- Assign positions: x = rapidity, y = layer index or parent depth ---
        # Compute "depth" from the number of ancestors
        depth = {}
        for node in G.nodes:
            d = 0
            current = node
            visited = set()  # prevent cycles
            while True:
                # convert node index to array index (assuming nodes go 1..n)
                array_idx = current - 1  
                current_parent = event_dict["parent"][array_idx]
                if current_parent == 0 or current_parent is None:
                    break
                if current in visited:
                    # cycle detected, stop to avoid infinite loop
                    print(f"Warning: cycle detected at node {current}")
                    break
                visited.add(current)
                d += 1
                current = current_parent
            depth[node] = d

        # Now create position dict
        pos = {}
        for i, idx in enumerate(event_dict["index"]):
            x = rapidities[i]
            y = -depth.get(idx, 0)  # vertical spacing by depth
            pos[idx] = (x, y)

    elif layout == "hierarchy":

        # Compute layers by parent relationships
        layers = {}
        for i, parent in zip(event_dict["index"], event_dict["parent"]):
            if parent == 0:
                layers[i] = 0  # root
            else:
                layers[i] = layers.get(parent, 0) + 1

        # Assign to graph
        nx.set_node_attributes(G, layers, "subset")

        # Hierarchical layout (for near-trees)
        pos = nx.multipartite_layout(G, subset_key="subset", scale=1.)
    else:
        # Spring layout handles general DAGs well
        pos = nx.spring_layout(G, k=0.8, iterations=100, seed=42)

    # ---- Define PID-based styles ----
    def get_style(pid):
        # Try classifying the particle with 'particle' package
        try:
            part = Particle.from_pdgid(pid)
            if len(part.quarks) == 2:
                return {'color': '#FF8C00', 'shape': 'D', 'category': 'Meson'}    # orange diamond
            if len(part.quarks) == 3 or part.is_unflavoured_meson:
                return {'color': '#FF6F61', 'shape': 'p', 'category': 'Baryon'}   # coral red pentagon
        except Exception:
            pass
        abs_pid = abs(pid)
        if abs_pid in [22, 23, 24, 25]:  # photons, W, Z, Higgs
            return {'color': 'gold', 'shape': 'o'}
        elif abs_pid in range(1, 7):  # quarks
            return {'color': 'lightcoral', 'shape': 's'}
        elif abs_pid in [11, 13, 15, 12, 14, 16]:  # leptons
            return {'color': 'skyblue', 'shape': '^'}
        elif abs_pid == 21:  # gluon
            return {'color': 'lightgreen', 'shape': 'h'}
        else:  # default / unknown
            return {'color': 'gray', 'shape': 'o'}

    # ---- Draw each group separately ----
    f, ax = plt.subplots(figsize=(16, 12))
    for shape in ['o', 's', '^', 'h', 'D', 'p']:
        nodes = [n for n, d in G.nodes(data=True) if get_style(d[bank_prefix+'pid'])['shape'] == shape]
        colors = [get_style(G.nodes[n][bank_prefix+'pid'])['color'] for n in nodes]
        labels = {}
        for n in nodes:
            pid = G.nodes[n][bank_prefix+'pid']
            try:
                part = Particle.from_pdgid(pid)
                label = rf"${part.latex_name}$"
            except Exception:
                label = f"{pid}"  # fallback if unknown PID
            labels[n] = label
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors, node_shape=shape, node_size=900)
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)

    nx.draw_networkx_edges(
        G, pos,
        arrows=True,
        arrowstyle='-|>',
        arrowsize=15,
        connectionstyle="arc3,rad=0.0",
        min_source_margin=20,
        min_target_margin=20
    )

    if layout == "rapidities":
        # --- CFR/TFR boundary line ---
        ax = plt.gca()
        ax.axvline(x=0, color='k', linestyle='--', linewidth=1.5)

        # --- Label axes and legend ---
        ax.set_xlabel("Rapidity $y$ ($\gamma^*N$ frame)",usetex=True)
        ax.set_ylabel("Decay depth",usetex=True)
        ax.set_title(f"Event {event_idx}: Decay Graph with CFR/TFR Rapidity Separation",usetex=True)


    # ---- Legend ----
    legend_elements = [
        mlines.Line2D([], [], color='none', marker='o', markerfacecolor='gold', label='Boson', markersize=12),
        mlines.Line2D([], [], color='none', marker='s', markerfacecolor='lightcoral', label='Quark', markersize=12),
        mlines.Line2D([], [], color='none', marker='^', markerfacecolor='skyblue', label='Lepton', markersize=12),
        mlines.Line2D([], [], color='none', marker='h', markerfacecolor='lightgreen', label='Gluon', markersize=12),
        mlines.Line2D([], [], color='none', marker='D', markerfacecolor='#FF8C00', markersize=12, label='Meson'),    # orange diamond
        mlines.Line2D([], [], color='none', marker='p', markerfacecolor='#FF6F61', markersize=12, label='Baryon'),   # coral red pentagon
        mlines.Line2D([], [], color='none', marker='o', markerfacecolor='gray', label='Other', markersize=12)
    ]
    if layout == "rapidities":
        legend_elements.append(
            mlines.Line2D([], [], color='k', linestyle='--', label='Rapidity boundary (y=0)')
        )
    plt.legend(handles=legend_elements, title="Particle Types", loc="lower left", frameon=False)

    if layout != "rapidities":
        plt.title("Decay Chain (PID-colored)",usetex=True)
        plt.axis('off')
    plt.tight_layout()
    f.savefig(f'event_{layout}_{event_idx}.pdf')

def plot_batch_decay_chains(batch_dict,bank_prefix="MC::Lund_",batch_idx=0,layouts=["hierarchy","rapidities","mothers_only","default"]):
    """
    Given a batch dictionary with shape (batch_size, nparticles[event]),
    loop over each event and plot its decay chain.
    
    Parameters
    ----------
    batch_dict : dict
        Keys: 'pid', 'index', 'parent'
        Each value: list of lists (batch_size × nparticles[event])
    """
    n_events = len(batch_dict['pid'])
    for i in range(n_events):
        print(f"\nPlotting event {i+1}/{n_events}...")

        # Extract this event’s data
        event_dict = {
            bank_prefix+'pid': batch_dict[bank_prefix+'pid'][i],
            bank_prefix+'index': batch_dict[bank_prefix+'index'][i],
            bank_prefix+'parent': batch_dict[bank_prefix+'parent'][i],
            bank_prefix+'daughter': batch_dict[bank_prefix+'daughter'][i],
            bank_prefix+'px': batch_dict[bank_prefix+'px'][i],
            bank_prefix+'py': batch_dict[bank_prefix+'py'][i],
            bank_prefix+'pz': batch_dict[bank_prefix+'pz'][i],
            bank_prefix+'E': batch_dict[bank_prefix+'E'][i],

        }

        # Boost event dictionary
        beam_lepton_idx = 0 #NOTE: This is a raw python index not the lund index.
        scattered_lepton_idx = 3 #NOTE: This is a raw python index not the lund index.
        beam_lepton = (
            event_dict["E"][beam_lepton_idx],
            event_dict["px"][beam_lepton_idx],
            event_dict["py"][beam_lepton_idx],
            event_dict["pz"][beam_lepton_idx],
        )
        scattered_lepton = (
            event_dict["E"][scattered_lepton_idx],
            event_dict["px"][scattered_lepton_idx],
            event_dict["py"][scattered_lepton_idx],
            event_dict["pz"][scattered_lepton_idx],
        )
        boosted_event_dict, _ = boost_event_dict_to_gammaN(
            event_dict, beam_lepton, scattered_lepton, target_mass=Particle.from_pdgid(2212).mass/1000.0
        )
        event_dict.update(boosted_event_dict)

        # Set the event index
        event_idx = n_events * batch_idx + i

        # Plot it
        if isinstance(layouts,list):
            for layout in layouts:
                plot_decay_chain(event_dict,bank_prefix=bank_prefix,event_idx=event_idx,layout=layout)
        else:
            plot_decay_chain(event_dict,bank_prefix=bank_prefix,event_idx=event_idx,layout=layouts)


# ---------------------------------------------------------------------------
# Main script
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Plot decay chains from HIPO files using hipopy.")
    parser.add_argument(
        "--files",
        nargs="+",
        required=True,
        help="List of input HIPO filenames or glob patterns (e.g. misc/tes*.hipo)."
    )
    parser.add_argument(
        "--bank",
        default="MC::Lund",
        help="Name of the bank containing event information (default: NEW::bank)."
    )
    parser.add_argument(
        "--step",
        type=int,
        default=10,
        help="Batch size for event iteration (default: 10)."
    )
    parser.add_argument(
        "--max_events",
        type=int,
        default=10,
        help="Max number of events to analyze (default: 10)."
    )
    parser.add_argument(
        "--layouts",
        nargs="+",
        required=False,
        choices = ["hierarchy","rapidities","mothers_only","default"],
        default = ["hierarchy","rapidities","mothers_only","default"],
        help="Plot layouts."
    )

    args = parser.parse_args()

    print("#----------------------------------------------------------------------#")
    print("Reading HIPO files:", args.files)
    print("Using bank:", args.bank)
    print("Batch step size:", args.step)

    banks = [args.bank]

    for batch_idx, batch in enumerate(hp.iterate(args.files, banks, step=args.step)):
        print(batch.keys())

        # Convert awkward arrays to normal Python lists
        if f"{args.bank}_pid" in batch and f"{args.bank}_index" in batch and \
        f"{args.bank}_parent" in batch and f"{args.bank}_daughter" in batch:
            batch_dict = {
                "pid": ak.to_list(batch[f"{args.bank}_pid"]),
                "index": ak.to_list(batch[f"{args.bank}_index"]),
                "parent": ak.to_list(batch[f"{args.bank}_parent"]),
                "daughter": ak.to_list(batch[f"{args.bank}_daughter"]),
                "px": ak.to_list(batch[f"{args.bank}_px"]),
                "py": ak.to_list(batch[f"{args.bank}_py"]),
                "pz": ak.to_list(batch[f"{args.bank}_pz"]),
                "E": ak.to_list(batch[f"{args.bank}_energy"]),
            }
            plot_batch_decay_chains(batch_dict,bank_prefix="",batch_idx=batch_idx,layouts=args.layouts)
        else:
            print(f"Warning: Missing required fields in {args.bank}")

        if (batch_idx+1) % args.step == 0:
            print("counter =", (batch_idx+1))
        if (batch_idx+1)*args.step >= args.max_events:
            break

# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import networkx as nx
from particle import Particle
import hipopy.hipopy as hp
import argparse
import awkward as ak

def plot_decay_chain(event_dict,bank_prefix="MC::Lund_",event_idx=0):
    """
    Plots the decay chain of a particle physics event with PID-based node styles.

    Parameters
    ----------
    event_dict : dict
        Must contain:
            - bank_prefix+'pid'    : list or array of particle IDs
            - bank_prefix+'index'  : list or array of 1-based indices (1..nparticles)
            - bank_prefix+'parent' : list or array of parent indices (1..nparticles)
                          or 0 / None / -1 for root particle.
    """
    pids = event_dict[bank_prefix+'pid']
    indices = event_dict[bank_prefix+'index']
    parents = event_dict[bank_prefix+'parent']

    G = nx.DiGraph()

    for idx, pid in zip(indices, pids):
        G.add_node(idx, pid=pid)

    for idx, parent_idx in zip(indices, parents):
        if parent_idx not in [None, 0, -1]:
            G.add_edge(parent_idx, idx)

    # ---- Hierarchical layout ----
    def hierarchy_pos(G, root=None, width=1., vert_gap=0.3, vert_loc=0, xcenter=0.5):
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

    # ---- Define PID-based styles ----
    def get_style(pid):
        # Try classifying the particle with 'particle' package
        try:
            part = Particle.from_pdgid(pid)
            if len(part.quarks) == 2:
                return {'color': '#FF8C00', 'shape': 'D', 'category': 'Meson'}    # orange diamond
            if len(part.quarks) == 3:
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
        nx.draw_networkx_nodes(G, subgraph_positions, nodelist=nodes, node_color=colors, node_shape=shape, node_size=900)
        nx.draw_networkx_labels(G, subgraph_positions, labels=labels, font_size=9)

    nx.draw_networkx_edges(
        G, subgraph_positions,
        arrows=True,
        arrowstyle='-|>',
        arrowsize=15,
        connectionstyle="arc3,rad=0.0",
        min_source_margin=20,
        min_target_margin=20
    )

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
    plt.legend(handles=legend_elements, title="Particle Types", loc="lower left", fontsize=10)

    plt.title("Decay Chain (PID-colored)")
    plt.axis('off')
    plt.tight_layout()
    f.savefig(f'event_{event_idx}.pdf')

def plot_batch_decay_chains(batch_dict,bank_prefix="MC::Lund_",batch_idx=0):
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
            bank_prefix+'parent': batch_dict[bank_prefix+'parent'][i]
        }

        # Set the event index
        event_idx = n_events * batch_idx + i

        # Plot it
        plot_decay_chain(event_dict,bank_prefix=bank_prefix,event_idx=event_idx)


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

    args = parser.parse_args()

    print("#----------------------------------------------------------------------#")
    print("Reading HIPO files:", args.files)
    print("Using bank:", args.bank)
    print("Batch step size:", args.step)

    banks = [args.bank]

    for batch_idx, batch in enumerate(hp.iterate(args.files, banks, step=args.step)):
        print(batch.keys())

        # Convert awkward arrays to normal Python lists
        if f"{args.bank}_pid" in batch and f"{args.bank}_index" in batch and f"{args.bank}_parent" in batch:
            batch_dict = {
                "pid": ak.to_list(batch[f"{args.bank}_pid"]),
                "index": ak.to_list(batch[f"{args.bank}_index"]),
                "parent": ak.to_list(batch[f"{args.bank}_parent"]),
            }
            plot_batch_decay_chains(batch_dict,bank_prefix="",batch_idx=batch_idx)
        else:
            print(f"Warning: Missing required fields in {args.bank}")

        if (batch_idx+1) % args.step == 0:
            print("counter =", (batch_idx+1))
        if (batch_idx+1)*args.step >= args.max_events:
            break

# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

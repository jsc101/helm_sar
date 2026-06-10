"""Token-flattening helpers — turn a HELMObject into per-position symbol/label/descriptor
lists for the alignment table. Main chain and sidechain monomers are kept distinct."""
from __future__ import annotations


def main_chain_syms(obj) -> list[str]:
    """Main chain symbols only — used for NW alignment (no sidechain monomers)."""
    if not hasattr(obj, 'get_jpv_flat'):
        return []
    return [t['symbol'] for t in obj.get_jpv_flat() if t['is_main']]


def main_chain_labels(obj) -> list[str]:
    """Position labels for main chain only, parallel to main_chain_syms."""
    if not hasattr(obj, 'get_jpv_flat'):
        return []
    return [t['label'] for t in obj.get_jpv_flat() if t['is_main']]


def token_descs(obj) -> tuple[list[dict], dict]:
    """
    Return (main_descs, sc_desc_map) where:
      main_descs   — per-position descriptor dict for each main chain token (in order)
      sc_desc_map  — {main_pos: [desc, desc, ...]} sidechain descriptor lists (outward order)
    """
    if not hasattr(obj, 'get_jpv_flat'):
        return [], {}

    chains = {c['chain_id']: c for c in obj.data.get('_chains', [])}
    primary = obj.get_chain()
    if primary is None:
        return [], {}
    pid = primary['chain_id']

    all_chain_descs = {cid: obj.all_monomer_descriptors(cid) for cid in chains}
    primary_descs = all_chain_descs.get(pid, {})

    # Build sidechain position map (main_pos → list of (sc_chain_id, sc_helm_pos) outward)
    conn_graph = obj.data.get('connectivity_graph', [])
    sc_order: dict[int, list] = {}
    for conn in conn_graph:
        fc, fp = conn['from_chain'], conn['from_pos']
        tc, tp = conn['to_chain'],   conn['to_pos']
        if fc == pid and tc != pid:
            sc_id, main_attach, sc_attach = tc, fp, tp
        elif tc == pid and fc != pid:
            sc_id, main_attach, sc_attach = fc, tp, fp
        else:
            continue
        sc = chains.get(sc_id)
        if sc is None:
            continue
        n = len(sc['monomers'])
        positions = list(range(1, n + 1))
        if sc_attach == n:
            positions = list(reversed(positions))
        elif sc_attach != 1:
            idx = sc_attach - 1
            positions = positions[idx:] + list(reversed(positions[:idx]))
        sc_order[main_attach] = [(sc_id, p) for p in positions]

    main_descs = []
    sc_desc_map: dict[int, list] = {}
    for token in obj.get_jpv_flat():
        mp = token['main_pos']
        if token['is_main']:
            d = primary_descs.get(mp, {})
            main_descs.append(d.get('descriptors', {}))
        else:
            k = token['sc_pos'] - 1
            entry_list = sc_order.get(mp, [])
            if k < len(entry_list):
                sc_cid, sc_p = entry_list[k]
                d = all_chain_descs.get(sc_cid, {}).get(sc_p, {})
                sc_desc_map.setdefault(mp, []).append(d.get('descriptors', {}))
            else:
                sc_desc_map.setdefault(mp, []).append({})

    return main_descs, sc_desc_map

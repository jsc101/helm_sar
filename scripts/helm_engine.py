from __future__ import annotations

import re
import json
from collections import deque


class HELMObject:
    """
    Parsed HELM molecule. Constructed by HELMParser — do not instantiate directly.

    Primary data lives in .data:
      _chains           list of {chain_id, type, monomers[{pos, symbol, entry}]}
      connectivity_graph list of connection dicts (from_chain/pos/rgroup, to_chain/pos/rgroup, bond_type)
      connectivity_map  list of (int, int, bond_type) global-index tuples
      topology_class    str
      helm_string       original input

    Each monomer entry in _chains carries the MonomerDB entry dict under 'entry'
    (None if the symbol was not found in the DB).
    """

    def __init__(self, helm_json: dict) -> None:
        self.data = helm_json
        # flat global connectivity for graph operations
        self.connectivity = helm_json.get('connectivity_map', [])

    # ------------------------------------------------------------------ #
    # Chain / monomer access
    # ------------------------------------------------------------------ #

    def get_chain(self, chain_id: str | None = None) -> dict | None:
        """Return chain dict. None → primary (longest PEPTIDE) chain."""
        chains = self.data.get('_chains', [])
        if chain_id is not None:
            return next((c for c in chains if c['chain_id'] == chain_id), None)
        peptide = [c for c in chains if c['type'] == 'PEPTIDE']
        pool = peptide if peptide else chains
        return max(pool, key=lambda c: len(c['monomers'])) if pool else None

    def get_monomer(self, pos: int, chain_id: str | None = None) -> dict | None:
        """Return monomer dict at 1-based position (includes 'entry' from MonomerDB)."""
        chain = self.get_chain(chain_id)
        if chain is None:
            return None
        return next((m for m in chain['monomers'] if m['pos'] == pos), None)

    def positions(self, chain_id: str | None = None) -> list[int]:
        """Return sorted list of 1-based positions in the chain."""
        chain = self.get_chain(chain_id)
        if chain is None:
            return []
        return [m['pos'] for m in chain['monomers']]

    # ------------------------------------------------------------------ #
    # RDKit — monomer-level descriptors
    # ------------------------------------------------------------------ #

    def to_smiles(self, chain_id: str | None = None, cap_c: str = 'amide') -> str | None:
        """Assemble a chain to a fully capped canonical SMILES (delegates to rdkit_bridge)."""
        from scripts.rdkit_bridge import helm_obj_to_smiles
        return helm_obj_to_smiles(self, chain_id=chain_id, cap_c=cap_c)

    def monomer_mol(self, pos: int, chain_id: str | None = None):
        """
        Return an RDKit Mol for the monomer at *pos* with R-groups capped.
        Requires rdkit_bridge.monomer_to_mol.
        """
        monomer = self.get_monomer(pos, chain_id)
        if monomer is None:
            return None
        entry = monomer.get('entry')
        if entry is None:
            return None
        from scripts.rdkit_bridge import monomer_to_mol
        return monomer_to_mol(entry)

    def monomer_descriptors(self, pos: int, chain_id: str | None = None) -> dict:
        """
        Calculate RDKit physicochemical descriptors for the monomer at *pos*.
        Returns {} if the monomer is unknown or RDKit fails.
        """
        mol = self.monomer_mol(pos, chain_id)
        if mol is None:
            return {}
        from scripts.rdkit_bridge import mol_to_descriptors
        return mol_to_descriptors(mol)

    def all_monomer_descriptors(self, chain_id: str | None = None) -> dict[int, dict]:
        """
        Descriptors for every position in the chain.

        Returns
        -------
        dict[pos → {symbol, descriptors}]
        """
        chain = self.get_chain(chain_id)
        if chain is None:
            return {}
        result = {}
        for m in chain['monomers']:
            result[m['pos']] = {
                'symbol': m['symbol'],
                'in_db': m.get('entry') is not None,
                'descriptors': self.monomer_descriptors(m['pos'], chain_id),
            }
        return result

    def position_descriptors(
        self,
        positions: list[int],
        chain_id: str | None = None,
    ) -> dict[int, dict]:
        """
        Descriptors for a specific list of positions.

        Example
        -------
        obj.position_descriptors([1, 5, 10])
        """
        return {
            pos: {
                'symbol': (m := self.get_monomer(pos, chain_id)) and m['symbol'] or '?',
                'in_db': bool(m and m.get('entry')),
                'descriptors': self.monomer_descriptors(pos, chain_id),
            }
            for pos in positions
        }

    # ------------------------------------------------------------------ #
    # Graph
    # ------------------------------------------------------------------ #

    def get_monomer_counts(self) -> list[int]:
        """BFS component sizes over the global connectivity graph."""
        chains = self.data.get('_chains', [])
        n = sum(len(c['monomers']) for c in chains)
        if n == 0:
            return []

        adj: dict[int, list[int]] = {i: [] for i in range(n)}
        for u, v, _ in self.connectivity:
            if v not in adj[u]:
                adj[u].append(v)
            if u not in adj[v]:
                adj[v].append(u)

        visited: set[int] = set()
        counts: list[int] = []
        for i in range(n):
            if i not in visited:
                count = 0
                queue = deque([i])
                visited.add(i)
                while queue:
                    curr = queue.popleft()
                    count += 1
                    for nb in adj[curr]:
                        if nb not in visited:
                            visited.add(nb)
                            queue.append(nb)
                counts.append(count)
        return counts

    # ------------------------------------------------------------------ #
    # JPV
    # ------------------------------------------------------------------ #

    def get_jpv(self) -> str:
        """
        JPV notation — unified algorithm, topology-class-agnostic.

        Delimiter: '-'
        Connections annotated as (bond_idx,rgroup_num) on each participating residue.
          bond_idx  — sequential integer: intra-chain bonds first, then inter-chain
          rgroup_num — R-group number at that residue's end (1=N-term, 2=C-term, 3=sidechain)
        Inter-chain, one attachment  → secondary chain injected inline as res(branch).
        Inter-chain, two+ attachments → bracketed linker: res([branch]) at first attachment.
        Main chain = longest PEPTIDE chain.
        """
        def _rg_num(rg: str) -> str:
            m = re.search(r'\d+', rg)
            return m.group() if m else '?'

        chains     = self.data.get('_chains', self.data.get('chains', []))
        conn_graph = self.data.get('connectivity_graph', [])

        if not chains:
            return ''

        peptide_chains = [c for c in chains if c['type'] == 'PEPTIDE']
        primary     = max(peptide_chains if peptide_chains else chains,
                          key=lambda c: len(c['monomers']))
        primary_id  = primary['chain_id']
        chain_by_id = {c['chain_id']: c for c in chains}

        intra = [c for c in conn_graph
                 if c['from_chain'] == primary_id and c['to_chain'] == primary_id]
        inter = [c for c in conn_graph
                 if not (c['from_chain'] == primary_id and c['to_chain'] == primary_id)]

        # Annotation map for primary chain residues: pos → [annotation_strings]
        primary_annot: dict[int, list[str]] = {}

        # Intra-chain bonds: numbered 1..len(intra)
        for bidx, conn in enumerate(intra, start=1):
            primary_annot.setdefault(conn['from_pos'], []).append(
                f'({bidx},{_rg_num(conn["from_rgroup"])})')
            primary_annot.setdefault(conn['to_pos'],   []).append(
                f'({bidx},{_rg_num(conn["to_rgroup"])})')

        # Inter-chain bonds: numbered len(intra)+1..
        # Group by secondary chain so we can annotate secondary residues and inject branches
        sec_attachments: dict[str, list[dict]] = {}
        for i, conn in enumerate(inter):
            bidx = len(intra) + i + 1
            fc, tc = conn['from_chain'], conn['to_chain']
            if fc == primary_id:
                sec_id   = tc
                pri_pos, pri_rg = conn['from_pos'], conn['from_rgroup']
                sec_pos, sec_rg = conn['to_pos'],   conn['to_rgroup']
            else:
                sec_id   = fc
                pri_pos, pri_rg = conn['to_pos'],   conn['to_rgroup']
                sec_pos, sec_rg = conn['from_pos'], conn['from_rgroup']

            primary_annot.setdefault(pri_pos, []).append(f'({bidx},{_rg_num(pri_rg)})')
            sec_attachments.setdefault(sec_id, []).append({
                'bidx': bidx, 'pri_pos': pri_pos,
                'sec_pos': sec_pos, 'sec_rg': sec_rg,
            })

        # Build each secondary chain's JPV (with its own connection annotations)
        branch_at: dict[int, list[str]] = {}
        for sec_id, atts in sec_attachments.items():
            sec_chain = chain_by_id.get(sec_id)
            if not sec_chain:
                continue

            sec_annot: dict[int, list[str]] = {}
            for att in atts:
                sec_annot.setdefault(att['sec_pos'], []).append(
                    f'({att["bidx"]},{_rg_num(att["sec_rg"])})')

            sec_tokens = [
                m['symbol'] + ''.join(sec_annot.get(m['pos'], []))
                for m in sec_chain['monomers']
            ]
            sec_jpv = '-'.join(sec_tokens)
            if len(atts) > 1:
                sec_jpv = f'[{sec_jpv}]'

            first_pri_pos = min(att['pri_pos'] for att in atts)
            branch_at.setdefault(first_pri_pos, []).append(sec_jpv)

        # Assemble primary chain
        result = []
        for m in primary['monomers']:
            pos = m['pos']
            tok = m['symbol'] + ''.join(primary_annot.get(pos, []))
            for branch in branch_at.get(pos, []):
                tok += f'({branch})'
            result.append(tok)

        return '-'.join(result)

    # ------------------------------------------------------------------ #
    # Sidechain / lipidation helpers
    # ------------------------------------------------------------------ #

    def get_lipidation_pos(self) -> int | None:
        """Return the 1-based main chain position bearing a secondary chain attachment."""
        primary = self.get_chain()
        if primary is None:
            return None
        pid = primary['chain_id']
        for conn in self.data.get('connectivity_graph', []):
            fc, fp = conn['from_chain'], conn['from_pos']
            tc, tp = conn['to_chain'],   conn['to_pos']
            if fc == pid and tc != pid:
                return fp
            if tc == pid and fc != pid:
                return tp
        return None

    def get_sidechain_string(self) -> str:
        """
        Dash-separated sidechain monomers reading from main chain outward.
        First token = monomer attached to backbone; last = distal end.
        Returns '' for single-chain molecules.
        """
        primary = self.get_chain()
        if primary is None:
            return ''
        pid = primary['chain_id']
        chains = {c['chain_id']: c for c in self.data.get('_chains', [])}
        for conn in self.data.get('connectivity_graph', []):
            fc, fp = conn['from_chain'], conn['from_pos']
            tc, tp = conn['to_chain'],   conn['to_pos']
            if fc == pid and tc != pid:
                sc = chains.get(tc)
                sc_attach = tp
            elif tc == pid and fc != pid:
                sc = chains.get(fc)
                sc_attach = fp
            else:
                continue
            if sc is None:
                continue
            syms = [m['symbol'] for m in sc['monomers']]
            n = len(syms)
            if sc_attach == 1:
                return '-'.join(syms)
            elif sc_attach == n:
                return '-'.join(reversed(syms))
            else:
                idx = sc_attach - 1
                return '-'.join(syms[idx:] + list(reversed(syms[:idx])))
        return ''

    def get_jpv_flat(self) -> list[dict]:
        """
        Flat token list: main chain with sidechain monomers inserted after their
        attachment position (reading outward from the backbone).

        Each token dict:
            symbol   : monomer symbol string
            label    : 'N' for main chain pos N; 'N.K' for K-th sidechain monomer at pos N
            is_main  : True for main chain residues
            main_pos : 1-based main chain position
            sc_pos   : sidechain index (0 for main chain tokens, 1+ for sidechain)
        """
        primary = self.get_chain()
        if primary is None:
            return []
        pid = primary['chain_id']
        chains = {c['chain_id']: c for c in self.data.get('_chains', [])}

        insertions: dict[int, list[str]] = {}
        for conn in self.data.get('connectivity_graph', []):
            fc, fp = conn['from_chain'], conn['from_pos']
            tc, tp = conn['to_chain'],   conn['to_pos']
            if fc == pid and tc != pid:
                sc = chains.get(tc)
                main_attach, sc_attach = fp, tp
            elif tc == pid and fc != pid:
                sc = chains.get(fc)
                main_attach, sc_attach = tp, fp
            else:
                continue
            if sc is None:
                continue
            syms = [m['symbol'] for m in sc['monomers']]
            n = len(syms)
            if sc_attach == 1:
                ordered = syms
            elif sc_attach == n:
                ordered = list(reversed(syms))
            else:
                idx = sc_attach - 1
                ordered = syms[idx:] + list(reversed(syms[:idx]))
            insertions[main_attach] = ordered

        result = []
        for m in primary['monomers']:
            pos = m['pos']
            result.append({'symbol': m['symbol'], 'label': str(pos),
                           'is_main': True, 'main_pos': pos, 'sc_pos': 0})
            for k, sym in enumerate(insertions.get(pos, []), start=1):
                result.append({'symbol': sym, 'label': f'{pos}.{k}',
                                'is_main': False, 'main_pos': pos, 'sc_pos': k})
        return result

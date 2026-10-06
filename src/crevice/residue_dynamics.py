"""Persistent physical contacts and aligned residue-motion associations.

Pair occupancy pools all contact types, so a changing chemical label does not
split a persistent physical contact. Results describe associations, not causality.
"""
from __future__ import annotations
import math
from collections import Counter

from .presentation import FigureText, annotate_option


class ResidueDynamics:
    """Accumulate residue contacts and aligned centroid motion over trajectory frames.

    Used by ``crevice trajectory`` (``residue_dynamics`` of
    :func:`crevice.trajectory.analyze_trajectory`) and ``cavity-trajectory``.
    Call :meth:`add` once per frame with that frame's residue network (all frames
    must have the same residues and contact definition), then :meth:`finish`.

    Parameters
    ----------
    lining_residues : iterable of str, default ()
        Residue IDs that line the reference region (user-supplied or from the
        reference boundary); they define ``reference_lining_contact``,
        ``contact_role`` and the nonlocal lining partners.
    min_occupancy : float, default 0.75
        Contacts present in fewer than this fraction of frames are not reported;
        must be in [0, 1].

    Raises
    ------
    ValueError
        If ``min_occupancy`` is outside [0, 1].
    """
    def __init__(self,*,lining_residues=(),min_occupancy=.75):
        if not math.isfinite(min_occupancy) or not 0<=min_occupancy<=1:
            raise ValueError('min_occupancy must be in [0,1]')
        self.lining=set(lining_residues);self.threshold=min_occupancy
        self.nodes=None;self.settings=None;self.contacts={};self.positions=[];self.times=[];self.indices=[]

    def add(self,frame,network,*,frame_index,time):
        """Add one frame's residue network and residue centroids.

        Each residue pair counts once per frame, using its shortest edge; the
        interaction labels of the frames are pooled, so a changing chemical label
        does not split a persistent physical contact.

        Parameters
        ----------
        frame : StructureFrame
            The (aligned) frame; kept for the call signature, the positions come from
            the network nodes.
        network : ResidueInteractionNetwork
            Residue network of this frame (node metadata ``position`` = residue
            centroid, Å).
        frame_index : int
            Source frame index, recorded in the report.
        time : float
            Frame time, ps.

        Raises
        ------
        ValueError
            If the residues or contact settings differ from the first frame, a
            lining residue is unknown, or a centroid is not finite.
        """
        import numpy as np
        nodes={n.id:n for n in network.nodes if n.kind=='residue'}
        settings={k:network.metadata.get(k) for k in ['cutoff','distance_metric','include_hydrogen','include_hetero','interaction_rules']}
        if self.nodes is None:
            if self.lining-set(nodes):raise ValueError('Reference lining contains unknown residue identities')
            self.nodes=nodes;self.settings=settings
        elif nodes.keys()!=self.nodes.keys() or settings!=self.settings:
            raise ValueError('Residue dynamics requires identical node identities and contact definitions in every frame')
        ids=sorted(nodes)
        xyz=np.asarray([nodes[k].metadata['position'] for k in ids])
        if not np.isfinite(xyz).all():raise ValueError('Residue positions must be finite')
        row=len(self.positions);self.positions.append(xyz);self.times.append(time);self.indices.append(frame_index)
        seen={}
        for edge in network.edges:
            if edge.source not in nodes or edge.target not in nodes:continue
            pair=tuple(sorted((edge.source,edge.target)))
            if pair not in seen or edge.distance<seen[pair].distance:seen[pair]=edge
        for pair,edge in seen.items():
            item=self.contacts.setdefault(pair,{'frames':[],'distance_sum':0.,'min_distance':float('inf'),'types':Counter(),'metadata':edge.metadata})
            item['frames'].append(row);item['distance_sum']+=edge.distance;item['min_distance']=min(item['min_distance'],edge.distance);item['types'][edge.interaction]+=1

    def finish(self,*,aligned,include_contact_frames=False):
        """Contact occupancies, centroid fluctuations and motion correlations.

        Parameters
        ----------
        aligned : bool
            Whether the frames were rigidly aligned; without alignment the RMSF and
            the dynamic cross-correlations (DCCM) are ``None``.
        include_contact_frames : bool, default False
            Add each contact's frame rows as ``contact_frame_rows`` (used by
            cavity-trajectory statistics).

        Returns
        -------
        dict
            ``status`` (``complete`` or ``no_network_frames``), frame count, indices
            and times, the contact definition, ``min_occupancy``, ``residues``
            (group, ``reference_lining``, ``centroid_rmsf_A``, nonlocal lining
            partners) and ``contacts`` (one row per pair at or above
            ``min_occupancy``: occupancy, frame counts, mean and minimum distance in
            Å, contact-type fractions, nonlocal and chain flags, ``contact_role``
            = ``lining_pair``, ``lining_partner`` or ``other``, whole-trajectory and
            half-trajectory DCCM of the residue centroids), plus definitions and an
            interpretation note. Correlations are descriptive, not causal.
        """
        import numpy as np
        if not self.positions:return {'status':'no_network_frames'}
        xyz=np.asarray(self.positions);n,r,_=xyz.shape;ids=sorted(self.nodes)
        fluctuations=xyz-xyz.mean(axis=0);variance=np.mean(np.sum(fluctuations**2,axis=2),axis=0)
        lookup={key:i for i,key in enumerate(ids)}
        def dccm(a,b,start=0,stop=None):
            if not aligned:return None
            x=xyz[start:stop,a];y=xyz[start:stop,b]
            if len(x)<3:return None
            x=x-x.mean(axis=0);y=y-y.mean(axis=0)
            vx=np.mean(np.sum(x*x,axis=1));vy=np.mean(np.sum(y*y,axis=1))
            if min(vx,vy)<1e-12:return None
            return float(np.clip(np.mean(np.sum(x*y,axis=1))/np.sqrt(vx*vy),-1,1))
        rows=[]
        for (left,right),item in sorted(self.contacts.items()):
            count=len(item['frames']);occupancy=count/n
            if occupancy<self.threshold:continue
            i,j=lookup[left],lookup[right]
            c=dccm(i,j);first=dccm(i,j,0,n//2);last=dccm(i,j,n//2,None)
            meta=item['metadata'];lining_contact=left in self.lining or right in self.lining
            rows.append({'source':left,'target':right,'occupancy':occupancy,'contact_frames':count,'eligible_frames':n,
                         'mean_distance_when_present_A':item['distance_sum']/count,'minimum_distance_A':item['min_distance'],
                         'contact_type_fractions':{k:v/count for k,v in sorted(item['types'].items())},
                         'nonlocal_contact':meta.get('nonlocal_contact'),
                         'sequence_separation_observed':meta.get('sequence_separation_observed'),
                         'different_chains':meta.get('different_chains'),
                         'reference_lining_contact':lining_contact,
                         'contact_role':'lining_pair' if left in self.lining and right in self.lining else 'lining_partner' if lining_contact else 'other',
                         'dccm':c,'dccm_first_half':first,'dccm_second_half':last,
                         'dccm_same_sign_in_halves':bool(first*last>0) if first is not None and last is not None else None})
        if include_contact_frames:
            for row in rows:
                row['contact_frame_rows']=list(self.contacts[(row['source'],row['target'])]['frames'])
        partners={key:[] for key in ids}
        for row in rows:
            if not row['nonlocal_contact']:continue
            for a,b in [(row['source'],row['target']),(row['target'],row['source'])]:
                if b in self.lining:partners[a].append({'residue':b,'occupancy':row['occupancy'],'dccm':row['dccm']})
        residue_rows=[{'residue':key,'group':self.nodes[key].metadata.get('group','unassigned'),
                       'reference_lining':key in self.lining,
                       'centroid_rmsf_A':float(np.sqrt(variance[i])) if aligned else None,
                       'nonlocal_lining_partners':partners[key]} for i,key in enumerate(ids)]
        return {'status':'complete','frame_count':n,'frame_indices':self.indices,'frame_times':self.times,
                'contact_definition':self.settings,'min_occupancy':self.threshold,'residues':residue_rows,'contacts':rows,
                'reference_lining_residues':sorted(self.lining),
                'alignment_applied':aligned,'motion_representation':'selected-atom residue geometric centroid (not center of mass)',
                'motion_uncertainty':'not estimated; whole-trajectory and split-half correlations are descriptive',
                'lining_definition':'user-supplied reference-boundary membership, not per-frame lining occupancy',
                'pair_occupancy_definition':'any contact type; denominator includes every frame with the identical residue set',
                'interpretation':'Persistent proximity and correlated motion identify candidates for follow-up, not functional control or causal information flow'}


@annotate_option
def write_dynamics_summary(report,prefix,*,dpi=300,annotate=None):
    """Write the persistent lining-partner contact table and figure.

    ``report`` (the ``residue_dynamics`` result of
    :func:`crevice.trajectory.analyze_trajectory`, written by ``crevice
    trajectory --network-json``) lists residue-residue contacts with their
    occupancy (the fraction of analysed frames in contact) and the aligned
    centroid dynamic cross-correlation (DCCM) over the whole trajectory and each
    half. The figure shows up to 20 nonlocal contacts between a
    reference-lining residue and a non-lining partner, ranked by occupancy and
    then by absolute correlation. Correlated motion is descriptive, not causal.

    Parameters
    ----------
    report : dict
        Residue-dynamics report with a ``contacts`` list.
    prefix : str or Path
        Output path stem; ``_contacts.csv``, ``_residues.csv`` and ``_partners.png`` are appended.
    dpi : int, default 300
        Figure resolution.
    annotate : bool, optional
        Draw the figure title and the second line of the correlation axis label.
        ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        ``contacts_csv``, ``residues_csv`` and ``partners_png`` paths.

    Outputs
    -------
    PREFIX_contacts.csv : table
        Every contact: source, target, occupancy, frames, mean distance (Å),
        roles and DCCM values.
    PREFIX_residues.csv : table
        One row per residue: ``residue``, ``group``, ``reference_lining``,
        ``centroid_rmsf_A`` (aligned centroid fluctuation, Å) and
        ``nonlocal_lining_partners`` (semicolon-separated residue IDs; their
        occupancy and correlation are in ``PREFIX_contacts.csv``).
    PREFIX_partners.png : figure
        "Fraction of analyzed frames in contact" bars and "Aligned motion
        correlation" (-1 to 1) for each contact pair (tick labels
        ``source ↔ target``).

    Colours and representation
    --------------------------
    Blue-grey (``#547c9b``) occupancy bars; magenta (``#a34477``) dots for the
    whole-trajectory correlation and translucent magenta bars spanning the
    first- and second-half correlations; grey zero line.
    """
    import csv
    from pathlib import Path
    from .figures import _pyplot,_save_figure
    prefix=Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    contacts=report.get('contacts',[])
    csv_path=prefix.with_name(prefix.name+'_contacts.csv')
    fields=['source','target','occupancy','contact_frames','eligible_frames','mean_distance_when_present_A',
            'nonlocal_contact','reference_lining_contact','contact_role','dccm','dccm_first_half','dccm_second_half']
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows({key:row.get(key) for key in fields} for row in contacts)
    residues_path=prefix.with_name(prefix.name+'_residues.csv')
    residue_fields=['residue','group','reference_lining','centroid_rmsf_A','nonlocal_lining_partners']
    with residues_path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=residue_fields);writer.writeheader()
        writer.writerows({**{key:row.get(key) for key in residue_fields[:4]},
                          'nonlocal_lining_partners':';'.join(v['residue'] if isinstance(v,dict) else str(v)
                                                              for v in row.get('nonlocal_lining_partners') or [])}
                         for row in report.get('residues',[]))
    selected=sorted((r for r in contacts if r.get('contact_role')=='lining_partner' and r['nonlocal_contact']),
                    key=lambda r:(-r['occupancy'],-abs(r['dccm'] or 0),r['source'],r['target']))[:20]
    text=FigureText(f"{len(selected)} persistent nonlocal contacts between reference-lining residues and non-lining partners: "
                    "contact occupancy and aligned motion correlation (dot whole trajectory, bar first to second half)")
    plt=_pyplot();fig,axes=plt.subplots(1,2,figsize=(12,7),sharey=True,layout='constrained')
    labels=[r['source']+' ↔ '+r['target'] for r in selected]
    axes[0].barh(labels,[r['occupancy'] for r in selected],color='#547c9b');axes[0].set_xlim(0,1.02)
    axes[0].set_xlabel('Fraction of analyzed frames in contact');axes[0].invert_yaxis()
    for i,row in enumerate(selected):
        if row['dccm'] is not None:
            axes[1].scatter(row['dccm'],i,color='#a34477',s=35)
        if row['dccm_first_half'] is not None and row['dccm_second_half'] is not None:
            axes[1].plot([row['dccm_first_half'],row['dccm_second_half']],[i,i],color='#a34477',alpha=.4,lw=3)
    axes[1].set_xlim(-1.03,1.03);axes[1].axvline(0,color='#aaaaaa',lw=.7)
    text.axis_label(axes[1].set_xlabel,'Aligned motion correlation','Dot: whole trajectory · line: half-trajectory values')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    if not selected:text.notice(axes[0],.5,.5,'No persistent nonlocal reference-lining contacts',transform=axes[0].transAxes,ha='center',wrap=True)
    text.suptitle(fig,'Persistent lining-partner contacts · descriptive motion associations',fontsize=13)
    png=prefix.with_name(prefix.name+'_partners.png');_save_figure(fig,png,dpi=dpi,text=text)
    return {'contacts_csv':str(csv_path),'residues_csv':str(residues_path),'partners_png':str(png)}

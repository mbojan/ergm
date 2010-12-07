#ifndef WTCHANGESTAT_H
#define WTCHANGESTAT_H

#include "wtedgetree.h"

typedef struct WtModelTermstruct {
  void (*d_func)(Edge, Vertex*, Vertex*, double*, struct WtModelTermstruct*, WtNetwork*);
  	void (*s_func)(struct WtModelTermstruct*, WtNetwork*);
        void (*t_func)(struct WtModelTermstruct*, WtNetwork*);
	double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
} WtModelTerm;


/****************************************************
 Macros to make life easier                         *
/* binomial coefficient macro: */
#define CHOOSE(n,r) ((n)<(r) ? (0) : (my_choose((double)(n),(int)(r)))) 

/* macros that tell whether a particular edge exists */
#define IS_OUTEDGE(a,b) (WtEdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define IS_INEDGE(a,b) (WtEdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define IS_UNDIRECTED_EDGE(a,b) IS_OUTEDGE(MIN(a,b), MAX(a,b))

/* macros that may be used to step through all in- or out-edges of a particular
   node.  These are used by the STEP_THROUGH_OUTEDGES and STEP_THROUGH_INEDGES 
   macros. */
#define MIN_OUTEDGE(a) (WtEdgetreeMinimum(nwp->outedges, (a)))
#define MIN_INEDGE(a) (WtEdgetreeMinimum(nwp->inedges, (a)))
#define NEXT_OUTEDGE(e) (WtEdgetreeSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE(e) (WtEdgetreeSuccessor(nwp->inedges,(e)))

/* macros to list each of the out-neighbors or in-neighbors, one at a time,
   of node a.  At each iteration of the loop, the variable v equals the node 
   number of the corresponding neighbor. */
#define STEP_THROUGH_OUTEDGES(a,e,v) for((e)=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES(a,e,v) for((e)=MIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE(e))
// Also execute for each edge, automatically adapting to undirected networks.
#define EXEC_THROUGH_OUTEDGES(a,e,v,subroutine) if(DIRECTED){ STEP_THROUGH_OUTEDGES(a,e,v) {subroutine} } else { EXEC_THROUGH_EDGES(a,e,v,subroutine) }
#define EXEC_THROUGH_INEDGES(a,e,v,subroutine) if(DIRECTED){ STEP_THROUGH_INEDGES(a,e,v) {subroutine} } else { EXEC_THROUGH_EDGES(a,e,v,subroutine) }
#define EXEC_THROUGH_EDGES(a,e,v,subroutine) { STEP_THROUGH_OUTEDGES(a,e,v) {subroutine}  STEP_THROUGH_INEDGES(a,e,v) {subroutine} }

// Non-adaptive versions. (I.e. ForceOUT/INEDGES.)
#define EXEC_THROUGH_FOUTEDGES(a,e,v,subroutine) STEP_THROUGH_OUTEDGES(a,e,v) {subroutine}
#define EXEC_THROUGH_FINEDGES(a,e,v,subroutine) STEP_THROUGH_INEDGES(a,e,v) {subroutine}

/* The OUTVAL and INVAL macros give the "other endnode" of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define OUTVAL(e) (nwp->outedges[(e)].value)
#define INVAL(e) (nwp->inedges[(e)].value)

/* macros for getting and setting the value of the (a,b) edge. */
#define GETWT(a,b) (WtGetEdge(a,b,nwp))
#define SETWT(a,b,w) (WtSetEdge(a,b,w,nwp))

#define N_NODES (nwp->nnodes) /* Total number of nodes in the network */
#define N_DYADS (nwp->directed_flag?(nnodes*(nnodes-1)):nwp->bipartite?nwp->bipartite*(nnodes-nwp->bipartite):(nnodes*(nnodes-1)/2))
#define OUT_DEG (nwp->outdegree) /* Vector of length N_NODES giving current outdegrees */
#define IN_DEG (nwp->indegree) /* Vector of length N_NODES giving current indegrees */
#define DIRECTED (nwp->directed_flag) /* 0 if network is undirected, 1 if directed */
#define N_EDGES (nwp->nedges) /* Total number of edges in the network currently */

/* 0 if network is not bipartite, otherwise number of first node of second type */
#define BIPARTITE (nwp->bipartite)

/* macros used for internal purposes:  assigning the next in- and out-edge when
   needed */
#define NEXT_INEDGE_NUM (nwp->next_inedge)
#define NEXT_OUTEDGE_NUM (nwp->next_outedge)

/* Vector of change statistics to be modified by the function*/
#define CHANGE_STAT (mtp->dstats)
/* Number of change statistics required by the current term */
#define N_CHANGE_STATS (mtp->nstats)

/* Vector of values passed via "inputs" from R */
#define INPUT_PARAM (mtp->inputparams)
#define N_INPUT_PARAMS (mtp->ninputparams) /* Number of inputs passed */

#define TOGGLEIND toggleind

/* macro to set all changestats to zero at start of function */
#define ZERO_ALL_CHANGESTATS() for(unsigned int TOGGLEIND=0; TOGGLEIND<N_CHANGE_STATS; TOGGLEIND++) CHANGE_STAT[TOGGLEIND]=0.0

/* macros to cycle through all toggles proposed for the current step, then
   make the current toggle in case of more than one proposed toggle, then
   undo all of the toggles to reset the original network state.  */
#define FOR_EACH_TOGGLE() for(unsigned int TOGGLEIND=0; TOGGLEIND<ntoggles; TOGGLEIND++)
// The idea here is to essentially swap the contents of the proposed
// weights with the current weights, and then swap them back when
// done.
#define HEAD head
#define TAIL tail
#define NEWWT newwt
#define OLDWT oldwt

#define GETOLDTOGGLEINFO() Vertex HEAD=heads[TOGGLEIND], TAIL=tails[TOGGLEIND]; double OLDWT=GETWT(HEAD,TAIL);
#define GETTOGGLEINFO() GETOLDTOGGLEINFO(); double NEWWT=weights[TOGGLEIND];

// SETWT_WITH_BACKUP(a) must be called _after_ GETTOGGLEINFO!
#define SETWT_WITH_BACKUP() {SETWT(HEAD,TAIL,NEWWT); weights[TOGGLEIND]=OLDWT;}
#define UNDO_SETWT() {GETOLDTOGGLEINFO(); SETWT(HEAD,TAIL,weights[TOGGLEIND]); weights[TOGGLEIND]=OLDWT;}
#define SETWT_IF_MORE_TO_COME() {if(TOGGLEIND+1<ntoggles) {SETWT_WITH_BACKUP()}}
#define UNDO_PREVIOUS_SETWTS() for(unsigned int TOGGLEIND=0; TOGGLEIND+1<ntoggles; TOGGLEIND++){UNDO_SETWT()}
// Brings together the above operations:
// For each toggle:
//    Get the current edge weight.
//    Calculate the change.
//    Back up the current edge weight by swapping weight[i] with current edge weight.
// For each toggle:
//    Undo the changes by swapping them back.
#define EXEC_THROUGH_TOGGLES(subroutine){ZERO_ALL_CHANGESTATS();FOR_EACH_TOGGLE(){ GETTOGGLEINFO(); {subroutine}; SETWT_IF_MORE_TO_COME();} UNDO_PREVIOUS_SETWTS();}

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define WtD_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *heads, Vertex *tails, double *weights, WtModelTerm *mtp, WtNetwork *nwp)
#define WtT_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtS_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)

WtD_CHANGESTAT_FN(d_from_s);
// This could be done more efficiently (saving a function call) 
// by assigning a function pointer as follows:
// #define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(*a)=d_from_s;
// However, it looks like it might confuse the function finding routines.
// In the future, it might be a good idea to have the initialization
// code autodetect when D_ function is not found, but S_ function is, and set it properly.
#define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(a){ d_from_s(ntoggles, heads, tails, weights, mtp, nwp); }

/* Not often used */
#define INPUT_ATTRIB (mtp->attrib)

#endif

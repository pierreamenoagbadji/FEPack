/*
  XLiFE++ is an extended library of finite elements written in C++
  Copyright (C) 2014  Lunéville, Eric; Kielbasiewicz, Nicolas; Lafranche, Yvon; Nguyen, Manh-Ha; Chambeyron, Colin

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*!
  \file Constraints.cpp
  \author E. Lunéville
  \since 20 jan 2014
  \date  20 jan 2014

  \brief Implementation of xlifepp::Constraints class functions
*/

#include "term.h"
#include "Constraints.hpp"

namespace xlifepp
{

// ---------------------------------------------------------------------------------------------------------------------------------
// constructor of Constraints object
// ---------------------------------------------------------------------------------------------------------------------------------
Constraints::Constraints(MatrixEntry* me, VectorEntry* ve)
  : matrix_p(me), rhs_p(ve), reduced(false), local(true), symmetric(true), isId(false) {}


// ---------------------------------------------------------------------------------------------------------------------------------
// //create ConstraintData from EssentialCondition, depending on the type of essential conditions
// ---------------------------------------------------------------------------------------------------------------------------------
Constraints::Constraints(const EssentialCondition& ec)
  : matrix_p(0), rhs_p(0), reduced(false), local(true), symmetric(true), isId(false)
{
  conditions_.push_back(ec);
  symmetric = false;

  if (ec.type()== _lfEc)   //essential condition given by a linear form
  {
    createLf(ec);      //from TermVector of a linearform
    return;
  }

  //special case u=0 with Argyris
  Space* sp = ec.unknown(1)->space();           //rootspace
  if(sp->interpolation()->type==_Argyris)
  {
      if(ec.diffOperator()==_id && ec.isHomogeneous()) {createArgyris(ec);return;}
      if(ec.diffOperator()==_ndotgrad && ec.isHomogeneous()) {createDirichlet(ec);return;}
      error("free_error","Argyris approximation supports only u=0 and ndotgrad(u)=0 boundary conditions");
  }

  // Lagrange u=0/op(g), Hdiv element u.n=0/op(g), Hrot element u^n=0/g^n
  if (ec.type()==_DirichletEc)
  {
    bool done = createDirichlet(ec);
    if (done) return;
  }

  //general method using interpolation on node of domain
  createNodal(ec);
}

// ---------------------------------------------------------------------------------------------------------------------------------
// copy constructor and assign operator
// ---------------------------------------------------------------------------------------------------------------------------------
Constraints::Constraints(const Constraints& c)
{
  matrix_p=0;
  rhs_p=0;
  copy(c);
}

Constraints& Constraints::operator=(const Constraints& c)
{
  if (this==&c) return *this;
  clear();
  copy(c);
  return *this;
}

void Constraints::copy(const Constraints& c)
{
  reduced=c.reduced;
  local=c.local;
  symmetric=c.symmetric;
  conditions_=c.conditions_;
  cdofsr_=c.cdofsr_;
  cdofsc_=c.cdofsc_;
  elcdofs_=c.elcdofs_;
  recdofs_=c.recdofs_;
  if (matrix_p!=0) delete matrix_p;
  if (rhs_p!=0) delete rhs_p;
  matrix_p=0; rhs_p=0;
  if (c.matrix_p!=0) matrix_p = new MatrixEntry(*c.matrix_p);
  if (c.rhs_p!=0) rhs_p = new VectorEntry(*c.rhs_p);
}

void Constraints::clear()
{
  if (matrix_p!=0) delete matrix_p;
  if (rhs_p!=0) delete rhs_p;
  matrix_p=0;
  rhs_p=0;
  cdofsr_.clear();
  cdofsc_.clear();
  elcdofs_.clear();
  recdofs_.clear();
  reduced=false;
}

// ---------------------------------------------------------------------------------------------------------------------------------
// //create ConstraintData from EssentialCondition, depending on the type of essential condition
// ---------------------------------------------------------------------------------------------------------------------------------
Constraints::~Constraints()
{
  if (matrix_p!=0) delete matrix_p;
  if (rhs_p!=0) delete rhs_p;
}

// ---------------------------------------------------------------------------------------------------------
/*! build Constraints set from EssentialConditions (list of essential conditions)
    - create Constraints object for each essential condition
    - merge constraints involving same unknown
    - if there exist a constraint coupling different unknowns merge all Constraints object in one Constraints object
    the output is a map of Constraints pointer indexed by unknown pointer;
    only one Constraints pointer in case of coupling conditions
*/
// ---------------------------------------------------------------------------------------------------------
std::map<const Unknown*, Constraints*> buildConstraints(const EssentialConditions& ecs)
{
  trace_p->push("buildConstraints");

  std::list<EssentialCondition>::const_iterator ite=ecs.begin();
  if (ite==ecs.end()) error("is_void", "ecs"); //void list

  //create one constraint for each essential condition
  std::vector<Constraints*> constraints(ecs.size());
  std::vector<Constraints*>::iterator itc=constraints.begin();
  for (; ite!=ecs.end(); ite++, itc++) *itc=new Constraints(*ite);

//  theCout<<"-------------------- build constraints -----------------------"<<eol;
//  for(itc=constraints.begin();itc!=constraints.end();++itc) theCout<<(**itc);

  //merge constraints
  std::map<const Unknown*, Constraints*>::iterator itm;
  std::map<const Unknown*, Constraints*> mconstraints =  mergeConstraints(constraints);

//  theCout<<"-------------------- merge constraints -----------------------"<<eol;
//  for(itm=mconstraints.begin();itm!=mconstraints.end();++itm) theCout<<(*itm->second);

  //reduce constraints
  for (itm=mconstraints.begin(); itm!=mconstraints.end(); itm++) itm->second->reduceConstraints();

//  theCout<<"-------------------- reduce constraints ----------------------"<<eol;
//  for(itm=mconstraints.begin();itm!=mconstraints.end();++itm) theCout<<(*itm->second);

  trace_p->pop();
  return mconstraints;
}

//=========================================================================================================
/*! construct constraints matrix in case of Dirichlet condition

    Lagrange scalar/vector unknown u :                      a*u = 0/g
    Raviart-Thomas (3D), Nedelec face (3D) scalar unknown : a*u.n = 0/g
    Nedelec (2D), Nedelec Edge (3D) scalar unknown :        a*u^n = 0
    Morlay scalar unknown :                                 a*u=g or grad(u).n=g

    build Id matrix and rhs vector 0/g
*/
//=========================================================================================================

bool Constraints::createDirichlet(const EssentialCondition& ec)
{
  //check if a standard Dirichlet condition
  const GeomDomain* dom=ec.domain(1);
  const Unknown* u=ec.unknown(1);
  Space* sp=u->space();           //rootspace
  if (!sp->isFE()) return false;
  FEType type =sp->interpolation()->type;
  bool direct = (type==_Lagrange      && ec.diffOperator()==_id)
                || (type==_RaviartThomas && ec.diffOperator()==_ndot && ec.isHomogeneous())
                || (type==_NedelecFace  && ec.diffOperator()==_ndot && ec.isHomogeneous())
                || (type==_NedelecEdge  && ec.diffOperator()==_ncross && ec.isHomogeneous())
                || ((type==_Morley || type==_Argyris)  && (ec.diffOperator()==_ndotgrad || ec.diffOperator()==_id));
  if (!direct) return false;  //go to general method

  trace_p->push("Constraints::createDirichlet");

  //get coefficient a
  complex_t a=ec.coefficient();    //coefficient of condition
  ValueType vta=_real;
  if (a.imag()!=0) vta=_complex;
  real_t ra=a.real();

  //create Id matrix
  Space* subsp=Space::findSubSpace(dom,sp);
  //if (subsp==0) subsp=new Space(*dom, *sp, sp->name()+"_"+dom->name());
  std::vector<number_t> dofids;
  if (subsp!=0) dofids = subsp->dofIds();              //dofs of subspace related to dom
  else // built dof using dofsOn function
  {
    const_cast<GeomDomain*>(dom)->updateParentOfSideElements();
    if(type==_Morley|| type==_Argyris) //special case of fourth order problem dealing with two types of essential conditions on a boundary, to be generalized in future
    {
        if(ec.diffOperator()==_id) dofids = sp->feSpace()->dofsOn(*dom,0,0,_noProjection); // value dof
        else if(ec.diffOperator()==_ndotgrad) dofids = sp->feSpace()->dofsOn(*dom,1,0,_dotnProjection);  // normal derivative dof
    }
    else dofids = sp->feSpace()->dofsOn(*dom);          //dofs located on dom
  }
  if (dofids.size()==0) error("dof_not_found");

  number_t n = dofids.size() * u->nbOfComponents();   //number of dofs
  matrix_p = new MatrixEntry(_idMatrix, _cs, _col, n, n, 1.);

  //create column dofs
  cdofsc_ = createCdofs(u,dofids);

  // create virtual constraint space and constraint unknown (row unknown)
  Function zero;                                                    //create void function to link it to a fake spectral space
  Space* csp=new Space(*dom, zero, n, 1,"C_"+unknownEcName(ec));    //create virtual constraint space
  const Unknown* csv= new Unknown(*csp, unknownEcName(ec), 1);      //create virtual constraint unknown
  cdofsr_=createCdofs(csv,csp->dofIds());                           //create constraint dof

  // create rhs vector entry
  const Function* fun=ec.funp();
  if (fun==0) rhs_p = new VectorEntry(real_t(0),n);
  else
  {
    ValueType vtfun=fun->valueType();
    rhs_p= new VectorEntry(vtfun,1,n);
    std::vector<Point> nodes(dofids.size());
    std::vector<Point>::iterator itp=nodes.begin();
    std::vector<number_t>::iterator itn=dofids.begin();
    for (; itn!=dofids.end(); ++itn, ++itp) *itp = sp->feSpace()->dofs[*itn-1].coords();
    std::vector<Vector<real_t> > ns;
    if (fun->normalRequired()) ns=computeNormalsAt(*dom,nodes); //compute normals on dom1
    if (vtfun==_real) buildRhs(fun,nodes,ns,real_t(0));
    else             buildRhs(fun,nodes,ns,complex_t(0));
    //divide by a if different from 1
    if (a != complex_t(1.,0.))
    {
      if (vtfun==_real && vta==_complex) {rhs_p->toComplex(); vtfun=_complex;}
      if (vtfun ==_complex) *rhs_p*=(1./a);
      else *rhs_p*=(1./ra);
    }

  }
  local=true;
  isId=true;
  symmetric=true;
  //end
  trace_p->pop();
  return true;
}

//=========================================================================================================
/*! General method to construct "nodal" constraint  a1*op1(u1)|d1 + a2*op2(u2)|d2 + ... = g
    where ai are complex scalars, opi are operators on unknowns ui and di are geometric domain
        - only one ui,di     :  u|d = g            Dirichlet condition
        - ui may be the same :  u|d1 - u|d2 = g    periodic condition
        - di may be the same :  u1|d - u2|d = g    transmission condition
        - ui,di are different:  u1|d1 - u2|d2 = g  generalized transmission condition
  When domains are different it is mandatory to give the maps mi relating d1 to di (see defineMap function),
  so the condition reads
                          a1*op1(u1)|d1 + a2*op2(u2)|m2(d1) + ... = g
  In any case, some nodes of d1 drive the computation :
                       a1*op1(u1)(Mk) + a2*op2(u2)(m2(Mk)) + ... = g(Mk)   Mk in d1
  The nodes Mk depends on the status of domain d1 :
    1) if domain d1 belongs to the mesh supporting u1 and u1 is a Lagrange interpolation
       then nodes used are all the nodes of the FE interpolation belonging to domain d1
    2) if domain d1 belongs to the mesh supporting u1 and u1 is a Nedelec or Raviart interpolation
       then nodes used are all the virtual coordinates of the dofs belonging to domain d1
    3) else nodes defining geometric domain are used

  Be cautious : u1|d1 + u2|d2 = g  is not equivalent to u2|d2 + u1|d1 = g

  The algorithm is the following
      - construct the list of nodes (Nk) of d1 regarding case
      - build matrix ai*opi(ui)|mi(di) for i=1,nbterms (see constraintsMatrix function)
      - concatenate them (see concatenateMatrix function)
      - construct right hand side

  \note the constraints matrix has to be a scalar matrix. It means that in case of vector unknown the constraints are split

                                                     | a  0  ...|
           Dirichlet and scalar unknown case     C = | 0 ...  0 |   diagonal matrix
                                                     |... 0   a |

                                                     | [a 0] [0 0]  ...  ...|     | a 0 0 0 ... ...|
                                                     | [0 a] [0 0]  ...  ...|     | 0 a 0 0 ... ...|
           Dirichlet and vector unknown(2D) case C = | [0 0] [a 0] [0 0] ...|  -> | 0 0 a 0 0 0 ...|
                                                     | [0 0] [0 a] [0 0] ...|     | 0 0 0 a 0 0 ...|
                                                     | [0 0] [0 0] [a 0] ...|     | 0 0 0 0 a 0 ...|
                                                     |  ...   ...   ...  ...|     | ... ... ... ...|
*/
//=========================================================================================================
void Constraints::createNodal(const EssentialCondition& ec)
{
  trace_p->push("Constraints::createNodal");
  if (ec.nbTerms()==0) error("is_void","term");

  const GeomDomain* dom1=ec.domain(1);
  const Unknown* u1=ec.unknown(1);
  if (dom1==0) error("null_pointer","domain");
  if (u1==0) error("null_pointer","unknown");
  Space* sp1=u1->space(), *subsp1=0;

  //built node list where apply constraints
  std::vector<number_t>::iterator itd;
  std::vector<Point>::iterator itp;
  std::vector<Point> nodes;
  FEType fet=sp1->feSpace()->interpolation()->type;
  if (dom1->mesh() == sp1->domain()->mesh() &&  fet ==_Lagrange)    //use nodes of u1 interpolation
  {
    subsp1=Space::findSubSpace(dom1,sp1);
    if (subsp1==0) subsp1=new Space(*dom1, *sp1, sp1->name()+"_"+dom1->name());
    if (ec.ecMethod!=_internalNodeEC)
    {
      std::vector<number_t> dofs=subsp1->dofIds();
      nodes.resize(dofs.size());
      itp=nodes.begin();
      for (itd=dofs.begin(); itd!=dofs.end(); ++itd, ++itp)
        *itp = subsp1->rootSpace()->feSpace()->dofs[*itd-1].coords();
    }
    else //use internal nodes : scaling of ref dof coordinates xi-> (1+n*xi)/(n+2) where n is the order of the element
    {
      const RefElement* refElt=subsp1->element(number_t(0)).refElt_p;
      number_t nbdofs=refElt->nbDofs(), nbelt=subsp1->nbOfElements(), p=0;
      number_t n=refElt->order();
      nodes.resize(nbelt*nbdofs); //assuming same element
      itp=nodes.begin();
      Point q0(std::vector<real_t>(refElt->dim(),1./real_t(n+2)));
      real_t alpha=real_t(n)/real_t(n+2);
      for (number_t k=0; k<nbelt; k++)
      {
        const Element& elt=subsp1->element(k);
        refElt = elt.refElt_p;
        GeomElement* gelt=elt.geomElt_p;
        MeshElement* melt = gelt->meshElement();
        if (melt==0) melt=gelt->buildSideMeshElement();
        if (melt->geomMapData_p == 0)  melt->geomMapData_p = new GeomMapData(melt);
        GeomMapData* mapdata=melt->geomMapData_p;
        const std::vector<RefDof*>& refDofs = refElt->refDofs;
        std::vector<RefDof*>::const_iterator itr=refDofs.begin();
        for (; itr!=refDofs.end(); ++itr, ++itp, p++)
        {
          Point q=(*itr)->point(); q*=alpha; q+=q0;
          *itp=mapdata->geomMap(q);
          number_t ns=nodes.size();
          if (p==nodes.size()) //increase the size of nodes
          {
            nodes.resize(ns+(nbelt-k)*nbdofs);
            itp=nodes.begin()+p;
          }
        }
      }
    }
  }
  else if (dom1->mesh() == sp1->domain()->mesh() && (fet==_NedelecEdge || fet==_NedelecFace || fet==_RaviartThomas || fet==_CrouzeixRaviart|| fet==_Morley))
  {
    //use internal side dofs "virtual location" given by dofs
    std::vector<number_t> dofnum=sp1->feSpace()->dofsOn(*dom1); //dofs located on dom1
    nodes.resize(dofnum.size());
    std::vector<Point>::iterator itp=nodes.begin();
    std::vector<number_t>::iterator itn=dofnum.begin();
    for (; itn!=dofnum.end(); ++itn, ++itp) *itp = sp1->feSpace()->dofs[*itn-1].coords();
  }
  else //default : use nodes from domain
  {
    std::set<number_t> nn=dom1->meshDomain()->nodeNumbers();
    const std::vector<Point>& meshnodes=dom1->mesh()->nodes;
    nodes.resize(nn.size());
    std::set<number_t>::iterator itn=nn.begin();
    itp=nodes.begin();
    for (; itn!=nn.end(); ++itn, ++itp)  *itp = meshnodes[*itn -1];
  }
  if (nodes.size()==0) error("is_void","nodes");
  if (subsp1==0) subsp1=sp1;

  //create matrix related to u1
  cit_opuval itc=ec.begin();
  matrix_p = constraintsMatrix(*itc->first, dom1, subsp1, itc->second, nodes, 0, cdofsc_);

  itc++;
  for (number_t i=1; i<ec.nbTerms(); ++i, ++itc)
  {
    std::vector<DofComponent> cdofsci;
    const Function* mi=0;
    const GeomDomain* domi=ec.domain(i+1);
    Space* spi=ec.unknown(i+1)->space(), *subspi=spi;
    SpaceType spti=spi->typeOfSpace();
    if (spti==_feSpace)
    {
        FEType feti=spi->feSpace()->interpolation()->type;
        if (domi->mesh() == spi->domain()->mesh() && (feti==_Lagrange))
        {
          subspi=Space::findSubSpace(domi,spi);
          if (subspi==0) subspi=new Space(*domi, *spi, spi->name()+"_"+domi->name());
        }
    }
    if (domi!=dom1)
    {
      mi=findMap(*dom1,*domi);
      if (mi == 0) error("domain_map_undefined",dom1->name(),domi->name());
    }
    MatrixEntry* mati = constraintsMatrix(*itc->first, domi, subspi, itc->second, nodes, mi, cdofsci);
    concatenateMatrix(*matrix_p, cdofsc_, *mati, cdofsci);
    delete mati;
  }

  // create virtual constraint space and constraint unknown (row unknown)
  Function zero;                                 //create void function to link it to a fake spectral space
  number_t nbr=matrix_p->nbOfRows();
  Space* csp=new Space(*dom1, zero, nbr, 1,"C_"+unknownEcName(ec)); //create virtual constraint space
  const Unknown* csv= new Unknown(*csp, unknownEcName(ec), 1);      //create virtual constraint unknown
  cdofsr_=createCdofs(csv,csp->dofIds());                           //create constraint dof

  // create rhs vector entry
  const Function* fun=ec.funp();
  if (fun==0) rhs_p = new VectorEntry(real_t(0),nbr);
  else
  {
    std::vector<Vector<real_t> > ns;
    if (fun->normalRequired()) ns=computeNormalsAt(*dom1,nodes); //compute normals on dom1
    ValueType vtfun=fun->valueType();
    //dimPair ds=fun->dims();             //in future, test compatibility of dimension with operator in constraints
    rhs_p= new VectorEntry(vtfun,1,nbr);
    if (vtfun==_real) buildRhs(fun,nodes,ns,real_t(0));
    else              buildRhs(fun,nodes,ns,complex_t(0));
  }

  local=(ec.nbTerms()==1 && ec.diffOperator()==_id);  //only u=g condition are considered as local condition
  trace_p->pop();
}

// ---------------------------------------------------------------------------------------------------------
/*! create constraints matrix from a term  a*op(u)|nodes where nodes is a collection of nodes :

         [ ... a*op(w_i1)(Mk)... a*op(w_i2)(Mk) .]  for node Mk

    the matrix is a scalar one and always stored in column sparse (csc)
    built also the column component dofs vector

    opu   : operator on unknown
    dom   : domain where opu is restricted
    spu   : space or subspace used for interpolation (may be different from space of unknown)
    coeff : coefficient applied to operator
    nodes : list of points
    fmap  : mapping function of nodes, may be 0 meaning id mapping

    cdofsc :  component dofs vector built by this routine

*/
// ---------------------------------------------------------------------------------------------------------
MatrixEntry* Constraints::constraintsMatrix(const OperatorOnUnknown& opu, const GeomDomain* dom, Space* spu, const complex_t& coeff,
    const std::vector<Point>& nodes, const Function* fmap,
    std::vector<DofComponent>& cdofsc)
{
  trace_p->push("Constraints::constraintsMatrix");

  const Unknown* u=opu.unknown();
  if (u==0) error("null_pointer","unknown");
  //if (u->space()->typeOfSpace() != _feSpace) error("not_fe_space_type", u->space()->name());

  //various initialization
  if (spu==0) spu=u->space();        //largest space if not defined
  dimen_t nbcu=u->nbOfComponents();
  dimen_t dimf=u->dimFun();
  if (nbcu > 1) dimf = nbcu;      // vector extension
  number_t n=nodes.size();
  dimen_t d,m;  //block size in operator computation
  const MeshDomain* mdom=spu->domain()->meshDomain();
  bool useParent = (!mdom->isSideDomain());
  dimen_t dimdom=dom->dim();

  //identify value type
  ValueType vt=_real;
  ValueType opuvt=opu.valueType();
  if (opuvt==_complex) vt=_complex;
  if (coeff.imag()!=0) vt=_complex;
  real_t rcoeff=coeff.real();

  SpaceType spt = spu->typeOfSpace();

  // special case of a spectral unknown, spectral basis sj
  // constraints matrix is a full matrix Cij = sj(Mi)
  // opu not yet managed !
  if (spt == _spSpace)
  {
    if (!opu.isId())
        error("free_error"," in Constraints::constraintsMatrix, spectral space case supports only Id operator");
    dimen_t dfsp = spu->dimFun();   // spectral fun dimension
    if (dfsp>1)
      error("free_error"," in Constraints::constraintsMatrix, vector spectral space not yet supported");
    number_t dimsp = spu->dimSpace();
    number_t nbcol = dimsp * dfsp; // number of cols
    ValueType vtsp = spu->valueType();
    if (vtsp==_complex && vt==_real) vt=_complex;
    std::vector<number_t> row0(n);
    std::vector<number_t>::iterator itn=row0.begin();
    for (number_t i=0;i<n; ++i, ++itn) *itn=i+1;
    std::vector<std::vector<number_t> > rows(nbcol,row0);
    MatrixStorage* ms= new ColCsStorage(n,nbcol,rows);
    MatrixEntry* mat= new MatrixEntry(vt,_scalar, ms);
    std::vector<Point>::const_iterator itp;
    number_t i=1;
    for (itp=nodes.begin(); itp!=nodes.end(); ++itp, ++i)
    {
      Point pt=*itp;
      if (fmap!=0) pt=(*fmap)(*itp, pt);  //map to
      // compute spectral function
      if (vtsp==_complex)
      {
        Vector<complex_t> spfs;
        spu->spSpace()->spectralFun(pt, spfs);
        Vector<complex_t>::iterator itf=spfs.begin();
        for (number_t k=1;k<=nbcol; ++k, ++itf ) mat->setEntry(i,k,coeff* *itf);
      }
      else // real spectral function
      {
        Vector<real_t> spfs;
        spu->spSpace()->spectralFun(pt, spfs);
        Vector<real_t>::iterator itf=spfs.begin();
        if (vt==_real)
        {
          for (number_t k=1;k<=nbcol; ++k, ++itf )
            mat->setEntry(i,k,rcoeff* *itf);
        }
        else
        {
          for (number_t k=1;k<=nbcol; ++k, ++itf )
            mat->setEntry(i,k,coeff* *itf);
        }
      }
    }
    std::vector<number_t> dofsp(dimsp);
    itn=dofsp.begin();
    for (number_t i=0;i<dimsp; ++i, ++itn) *itn=i+1;
    cdofsc = createCdofs(u,dofsp);
    trace_p->pop();
    return mat;
  }


  //create the map geomelement number -> element*
  std::map<number_t, const Element*> gelt2elt;
  for (number_t k=0; k<spu->nbOfElements(); k++)
  {
    const Element* elt=spu->element_p(k);
    gelt2elt[elt->geomElt_p->number()]=elt;
  }

  //loop initialization
  std::vector<VectorEntry> opuws(n);              //list of values for each node
  std::vector<std::vector<number_t> > dofs(n);    //list of dofs for each node
  std::map<number_t, number_t> dofmap;            //map global dof numbering to local dof numbering
  std::map<number_t, number_t>::iterator itm;
  std::vector<VectorEntry>::iterator itwu=opuws.begin();
  std::vector<std::vector<number_t> >::iterator itdu=dofs.begin();
  std::vector<Point>::const_iterator itp;
  Vector<real_t>* np=0;   //pointer to normal or other vector involved in operator
  number_t k=1, i=1;

  //main construction loop
  for (itp=nodes.begin(); itp!=nodes.end(); ++itp, ++itwu, ++itdu, ++i)
  {
    //geometric stuff
    Point pt=*itp;
    if (fmap!=0) pt=(*fmap)(*itp, pt);  //map to
    Point q=pt;
    real_t md; //distance to nearest element
    GeomElement* belt = dom->meshDomain()->nearest(pt,md);
    if (belt==0) error("point_not_on_boundary");
    if (dimdom > 0 && md > belt->measure()) error("point_too_far_from_boundary",md);
    const GeomElement* gelt = belt;
    if (useParent)
    {
      gelt = belt->parentInDomain(mdom);       //first parent in mdom
      if (gelt==0) gelt = mdom->nearest(q,md);  //try to locate element in mdom
    }
    const Element* elt=gelt2elt[gelt->number()];
    if (opu.normalRequired())
    {
      MeshElement* melt = belt->meshElement();
      if (melt==0) melt=belt->buildSideMeshElement();
      if (melt->geomMapData_p == 0)  melt->geomMapData_p = new GeomMapData(melt);
      GeomMapData* mapdata=melt->geomMapData_p;
      Point q=mapdata->geomMapInverse(pt);
      mapdata->computeOrientedNormal();
      np=&mapdata->normalVector;
    }

    //numbering stuff
    std::vector<number_t> dofn=elt->dofNumbers;
    std::vector<number_t>::iterator itd=dofn.begin(), itdc=itd;
    bool der1=opu.diffOrder()>0, der2=opu.diffOrder()>1;
    //compute shape functions and operator onto
    ShapeValues shv=elt->computeShapeValues(pt,der1,der2);
    if (nbcu > 1) shv.extendToVector(nbcu);                            //extend scalar shape functions to nbc vector shape functions
    if (vt==_real)  //all is real
    {
      Vector<real_t> val;
      if (opu.hasFunction()) opu.eval(pt, shv.w, shv.dw, shv.d2w, dimf, val, d, m, np);  //evaluate differential operator with function
      else                   opu.eval(shv.w, shv.dw, shv.d2w, dimf, val, d, m, np);      //evaluate differential operator
      //clean val
      Vector<real_t>::iterator itv=val.begin(), itvc=val.begin(), itv2;
      number_t nt=number_t(d*nbcu);   //number of values linked to a dof
      number_t nbd=0;
      for (itd=dofn.begin(); itd!=dofn.end(); ++itd)
      {
        real_t s=0.;
        itv2=itv;
        for (number_t j=0; j<nt; j++, ++itv) s+=std::abs(*itv);
        if (s>theTolerance) //keep values
        {
          for (number_t j=0; j<nt; j++, ++itvc, ++itv2) *itvc=*itv2;
          *itdc=*itd;
          itdc++;
          nbd++;
        }
      }
      dofn.resize(nbd);
      val.resize(nbd*nt);
      *itwu=VectorEntry(val); //store val as a VectorEntry
    }
    else //some are complex, store in complex
    {
      Vector<complex_t> val;
      if (opuvt==_real)
      {
        Vector<real_t> valr;
        if (opu.hasFunction()) opu.eval(pt, shv.w, shv.dw, shv.d2w, dimf, valr, d, m, np);  //evaluate differential operator with function
        else                  opu.eval(shv.w, shv.dw, shv.d2w, dimf, valr, d, m, np);      //evaluate differential operator
        val=cmplx(valr);
      }
      else
      {
        if (opu.hasFunction()) opu.eval(pt, shv.w, shv.dw, shv.d2w, dimf, val, d, m, np);  //evaluate differential operator with function
        else                  opu.eval(shv.w, shv.dw, shv.d2w, dimf, val, d, m, np);      //evaluate differential operator

      }
      //clean val
      Vector<complex_t>::iterator itv=val.begin(), itvc=itv, itv2;
      number_t nt=d*nbcu;   //number of values linked to a dof
      number_t nbd=0;
      for (itd=dofn.begin(); itd!=dofn.end(); ++itd)
      {
        real_t s=0.;
        itv2=itv;
        for (number_t j=0; j<nt; j++, ++itv) s+=std::abs(*itv);
        if (s>theTolerance) //keep values
        {
          for (number_t j=0; j<nt; j++, ++itvc, ++itv2) *itvc=*itv2;
          *itdc++=*itd;
          nbd++;
        }
      }
      dofn.resize(nbd);
      val.resize(nbd*nt);
      *itwu=VectorEntry(val);
    }

    //update dof numbering
    for (itd=dofn.begin(); itd!=dofn.end(); ++itd)
      if (dofmap.find(*itd)==dofmap.end()) dofmap[*itd] = k++;
    *itdu=dofn;

  } //end of loop

  //create cdofs vector
  std::vector<number_t> dofv(dofmap.size());
  std::vector<number_t>::iterator itd=dofv.begin(); k=1;
  for (itm=dofmap.begin(); itm!=dofmap.end(); ++itm, ++itd,++k)
  {
    itm->second=k;
    *itd=itm->first;
  }
  cdofsc = createCdofs(u,dofv);
  number_t nbcol = cdofsc.size();   //differs from nbdof if vector unknown

  //create row index (ordering using dofmap) and matrix storage
  std::vector<std::vector<number_t> > rows(nbcol);
  k=1;
  if (nbcu==1 && d==1)  //scalar case
  {
    for (itdu=dofs.begin(); itdu!=dofs.end(); ++itdu, ++k)
      for (itd=itdu->begin(); itd!=itdu->end(); ++itd) rows[dofmap[*itd]-1].push_back(k);
  }
  else //vector case
  {
    for (itdu=dofs.begin(); itdu!=dofs.end(); ++itdu, k+=d)    //loop on dofs
    {
      for (itd=itdu->begin(); itd!=itdu->end(); ++itd)
      {
        //number_t c = d*(dofmap[*itd]-1);
        number_t c = nbcu*(dofmap[*itd]-1);
        for (number_t i=0; i<d; ++i)
          for (number_t j=0; j<nbcu; ++j) rows[c+j].push_back(k+i);
      }
    }
  }

  //create storage and matrix
  number_t nbrow = n*d;
  MatrixStorage* ms= new ColCsStorage(nbrow,nbcol,rows);
  MatrixEntry* mat= new MatrixEntry(vt,_scalar, ms);

  //fill matrix (could be optimized)
  itwu=opuws.begin(); k=1;
  if (nbcu==1 && d==1)  //scalar case
  {
    for (itdu=dofs.begin(); itdu!=dofs.end(); ++itdu, ++itwu, ++k)
    {
      number_t i=1;
      for (itd=itdu->begin(); itd!=itdu->end(); ++itd, ++i)
      {
        if (vt==_real)
        {
          real_t r;
          itwu->getEntry(i,r);
          mat->setEntry(k,dofmap[*itd],rcoeff*r);
        }
        else
        {
          complex_t c;
          itwu->getEntry(i,c);
          mat->setEntry(k,dofmap[*itd],coeff*c);
        }
      }
    }
  }
  else //vector case
  {
    for (itdu=dofs.begin(); itdu!=dofs.end(); ++itdu, ++itwu, k+=d)
    {
      number_t l=1;
      for (itd=itdu->begin(); itd!=itdu->end(); ++itd)
      {
        for (number_t i=0; i<d; ++i)
        {
          for (number_t j=0; j<nbcu; ++j, ++l)
          {
            number_t r = k+i;
            number_t c = nbcu*(dofmap[*itd]-1)+j +1;
            if (vt==_real)
            {
              real_t vr;
              itwu->getEntry(l,vr);
              mat->setEntry(r,c,rcoeff*vr);
            }
            else
            {
              complex_t vc;
              itwu->getEntry(l,vc);
              mat->setEntry(r,c,coeff*vc);
            }
          }
        }
      }
    }
  }

  trace_p->pop();
  return mat;
}

// ---------------------------------------------------------------------------------------------------------------------------------
/*! create constraint data for a lf condition lf(u)=c (one row constraint), lf represented by a TermVector v
    the constraints matrix has to be a scalar matrix. It means that in case of vector unknown the constraints are split

           scalar unknown case     C = | v1 ...  vn |   row matrix

           vector unknown(2D) case C = | [v11 v12] [v21 v22] ...|  -> | v11 v12 v21 v22 ....|
*/
// ---------------------------------------------------------------------------------------------------------------------------------
void Constraints::createLf(const EssentialCondition& ec)
{
  trace_p->push("Constraints::createLf");
  if (ec.type()!=_lfEc) error("ec_bad_ectype",words("essential condition",ec.type()), words("essential condition",_lfEc));
  if (ec.lfp()==0) error("null_pointer","termVector (lf)");

  TermVector tv(*ec.lfp());
  tv.toGlobal();  //move to global representation
  VectorEntry* v=tv.actual_entries();
  number_t n=v->size();
  std::vector<number_t> un(1,1);
  std::vector<std::vector<number_t> > rowindices(n,un);
  MatrixStorage* ms= new ColCsStorage(1, n, rowindices);
  ValueType vt=v->valueType_;
  matrix_p = new MatrixEntry(vt,_scalar,ms);
  //copy entries
  if (vt==_real)
  {
    Vector<real_t>::const_iterator itv=v->rEntries_p->begin();
    Vector<real_t>::iterator itm=matrix_p->rEntries_p->values().begin();
    itm++;//jump first value
    for (; itv!=v->rEntries_p->end(); itv++, itm++) *itm=*itv;
  }
  else
  {
    Vector<complex_t>::iterator itv=v->cEntries_p->begin();
    Vector<complex_t>::iterator itm=matrix_p->cEntries_p->values().begin();
    itm++; //jump first value
    for (; itv!=v->cEntries_p->end(); itv++, itm++) *itm=*itv;
  }
  complex_t c=ec.clf;
  if (std::abs(c.imag())==0.)
    rhs_p = new VectorEntry(c.real(),1);
  else
    rhs_p = new VectorEntry(c,1);

  //create virtual constraint space
  const GeomDomain* dom=tv.firstSut()->domain(); //take arbitrary domain
  Function zero;
  Space* csp=new Space(*dom, zero, 1, 1,"C_"+unknownEcName(ec));
  const Unknown* csv= new Unknown(*csp, unknownEcName(ec), 1);
  //update cdofs
  cdofsc_=tv.cdofs();
  cdofsr_=createCdofs(csv,csp->dofIds());
  local=false;
  trace_p->pop();
}

// ---------------------------------------------------------------------------------------------------------
/*! create constraints for condition u=0 on gamma in case of Argyris approximations
    u=0 on gamma  -> u_i=0 but also null tangential derivatives at vertices :
    implies dy.u*nx-dx.u*ny = 0 -> dxu_i*nix-dyu_i*niy = 0
            dyy.u*nx*nx+dxx.u*ny*ny-2dxy.u*nx*ny =0 -> dyyu_i*nix*nix+dxxu_i*niy*niy-2dxyu_i*nix*niy =0
    with u_i,dxu_i,dyu_i, ... dof on vertex i
    NOTE : because side element are not defined for Argyris, we use more complex algorithm that travels geometric side elements
           to retry dof located on side gamma and then create constraints on value and derivative dofs
           because tangential derivatives at vertex may differ from side elements sharing vertex, both tangential derivative constraints are added
           and the reduction process will average them. This average process may induce non local constraints and thus matrix storage change
*/
// ---------------------------------------------------------------------------------------------------------
void Constraints::createArgyris(const EssentialCondition& ec)
{
  trace_p->push("Constraints::createNodal");
  if (ec.nbTerms()==0) error("is_void","term");
  local=false;

  const GeomDomain* dom=ec.domain(1);
  const Unknown* u=ec.unknown(1);
  if (dom==0) error("null_pointer","domain");
  if (u==0) error("null_pointer","unknown");
  Space* spu=u->space();
  std::map<GeomElement*, number_t>& gelt2elt = spu->feSpace()->gelt2elt;
  if(gelt2elt.size()==0) spu->buildgelt2elt();   //build map GeomElement to Element
  const std::vector<GeomElement*>& geoElts = dom->meshDomain()->geomElements;
  number_t nbelt=geoElts.size();
  //create 3*2*nbelt constraint in a sparse matrix (some are conflicting, will be reduced later)
  number_t nbrow=6*nbelt;
  SparseMatrix<real_t> sp(nbrow,1);
  std::vector<GeomElement*>::const_iterator itge;
  number_t k=1;
  for(itge = geoElts.begin(); itge != geoElts.end(); itge++)  //loop on geometric side element
  {
     std::vector<GeoNumPair>& parsides = (*itge)->parentSides(); //parents and side numbers
     std::vector<GeoNumPair>::iterator itvg=parsides.begin();
     for(; itvg!=parsides.end(); itvg++)   //travel parents
     {
       std::map<GeomElement*, number_t>::iterator itm=gelt2elt.find(itvg->first);
       if(itm!=gelt2elt.end())  //element found
       {
         const Element& elt = spu->feSpace()->elements[itm->second];
         const std::vector<number_t>& ds=elt.dofNumbers;
         number_t s = itvg->second;
         std::vector<real_t> n=itvg->first->normalVector(s);
         //theCout<<"side ="<<s<<", n="<<n<<", ds="<<ds<<eol;
         //  s=1 : dof 1-6 7-12   -> w1=w7=0   w2.ny-w3.nx=0   w8.ny-w9.nx=0     w6*nx*nx+w4*ny*ny-2w5*nx*ny =0  w12*nx*nx+w10*ny*ny-2w11*nx*ny =0
         //  s=2 : dof 7-12 13-18 -> w7=w13=0  w8.ny-w9.nx=0 w14.ny-w15.nx=0  w12*nx*nx+w10*ny*ny-2w11*nx*ny =0  w18*nx*nx+w16*ny*ny-2w17*nx*ny =0
         //  s=3 : dof 13-18 1-6  -> w1=w13=0  w2.ny-w3.nx=0 w14.ny-w15.nx=0     w6*nx*nx+w4*ny*ny-2w5*nx*ny =0  w18*nx*nx+w16*ny*ny-2w17*nx*ny =0
         if(s==1 || s==3)
         {
           sp(k,ds[0])=1.; k++;
           if(std::abs(n[1])>theEpsilon || std::abs(n[0])>theEpsilon )
           {
             if(std::abs(n[0])>theEpsilon)
             {
               sp(k,ds[2])=-n[0];
               sp(k+1,ds[5])= n[0]*n[0];
             }
             if(std::abs(n[1])>theEpsilon)
             {
               sp(k,ds[1])  = n[1];
               sp(k+1,ds[3])= n[1]*n[1];
             }
             if(std::abs(n[0]*n[1])>theEpsilon) sp(k+1,ds[4])=-2*n[0]*n[1];
             k+=2;
           }
         }
         if(s==1 || s==2)
         {
           sp(k,ds[6])=1.; k++;
           if(std::abs(n[1])>theEpsilon || std::abs(n[0])>theEpsilon )
           {
             if(std::abs(n[0])>theEpsilon)
             {
               sp(k,ds[8])=-n[0];
               sp(k+1,ds[11])= n[0]*n[0];
             }
             if(std::abs(n[1])>theEpsilon)
             {
               sp(k,ds[7])  = n[1];
               sp(k+1,ds[9])= n[1]*n[1];
             }
             if(std::abs(n[0]*n[1])>theEpsilon) sp(k+1,ds[10])=-2*n[0]*n[1];
             k+=2;
           }
         }
         if(s==2 || s==3)
         {
            sp(k,ds[12])=1.; k++;
            if(std::abs(n[1])>theEpsilon || std::abs(n[0])>theEpsilon )
           {
             if(std::abs(n[0])>theEpsilon)
             {
               sp(k,ds[14])=-n[0];
               sp(k+1,ds[17])= n[0]*n[0];
             }
             if(std::abs(n[1])>theEpsilon)
             {
               sp(k,ds[13])  = n[1];
               sp(k+1,ds[15])= n[1]*n[1];
             }
             if(std::abs(n[0]*n[1])>theEpsilon) sp(k+1,ds[16])=-2*n[0]*n[1];
             k+=2;
           }
         }
       }
     }
   }
   //dofs involved
   std::map<number_t,number_t> dofmap;
   k=1;
   std::map<NumPair,real_t>::iterator its;
   for(its=sp.begin();its!=sp.end();++its)
     if(dofmap.find(its->first.second)==dofmap.end()) dofmap[its->first.second]=k++;
   //theCout<<"sp="<<eol<<sp<<eol<<"dense sp="<<eol<<sp.toMatrix()<<eol<<"dofmap : "<<dofmap<<eol<<std::flush;
   number_t nbcol=dofmap.size();
   //create storage and matrix
   std::vector<std::vector<number_t> > rows(nbcol);
   for(its=sp.begin();its!=sp.end();++its)
       rows[dofmap[its->first.second]-1].push_back(its->first.first);
   //theCout<<"rows : "<<rows<<eol;
   MatrixStorage* ms= new ColCsStorage(nbrow,nbcol,rows);
   matrix_p = new MatrixEntry(_real,_scalar, ms);
   for(its=sp.begin();its!=sp.end();++its)
      matrix_p ->setEntry(its->first.first,dofmap[its->first.second],its->second);
   //create cdofsc vector
  std::vector<number_t> dofv(nbcol);
  std::vector<number_t>::iterator itd=dofv.begin();
  //for (std::map<number_t,number_t>::iterator itm=dofmap.begin(); itm!=dofmap.end(); ++itm, ++itd) *itd=itm->first;
  for (std::map<number_t,number_t>::iterator itm=dofmap.begin(); itm!=dofmap.end(); ++itm) dofv[itm->second-1]=itm->first;
  cdofsc_ = createCdofs(u,dofv);
  //theCout<<"cdofsc : "<<cdofsc_<<eol;
  // create virtual constraint space and constraint unknown (row unknown)
  Function zero;                                                      //create void function to link it to a fake spectral space
  Space* csp=new Space(*dom, zero, nbrow, 1,"C_"+unknownEcName(ec));  //create virtual constraint space
  const Unknown* csv= new Unknown(*csp, unknownEcName(ec), 1);        //create virtual constraint unknown
  cdofsr_=createCdofs(csv,csp->dofIds());                             //create constraint dof
  // create rhs vector entry (homogeneous condition)
  rhs_p = new VectorEntry(real_t(0),nbrow);

  trace_p->pop();
}

// ---------------------------------------------------------------------------------------------------------
/*! concatenate 2 constraint matrices mat1 and mat2 with no common column dofs
    mat1 and mat2 are scalar matrix stored in cs column, mat1 and cdofsc1 are updated
*/
void Constraints::concatenateMatrix(MatrixEntry& mat1, std::vector<DofComponent>& cdofsc1,
                                    const MatrixEntry& mat2, const std::vector<DofComponent>& cdofsc2)
{
  trace_p->push("Constraints::concatenateMatrix");

  if (mat1.storageType()!=_cs) error("storage_unexpected",words("storage type",_cs),words("storage type",mat1.storageType()));
  if (mat2.storageType()!=_cs) error("storage_unexpected",words("storage type",_cs),words("storage type",mat2.storageType()));
  if (mat1.accessType()!=_col) error("access_unexpected",words("access type",_col),words("access type",mat1.accessType()));
  if (mat2.accessType()!=_col) error("access_unexpected",words("access type",_col),words("access type",mat1.accessType()));

  ValueType vt1=mat1.valueType_, vt2=mat2.valueType_, vt=vt1;
  if (vt==_real && vt2==_complex) vt=_complex;
  if (vt==_complex && vt1==_real) mat1.toComplex();

  //concatenate storage
  ColCsStorage* sto1 = reinterpret_cast<ColCsStorage*>(mat1.storagep()),
                * sto2 = reinterpret_cast<ColCsStorage*>(mat2.storagep());
  number_t nnz1=sto1->size(), nnz2=sto2->size();
  number_t nc1=sto1->nbOfColumns(), nc2=sto2->nbOfColumns(), nc=nc1+nc2;

  std::vector<number_t>& colp1=sto1->colPointer();
  const std::vector<number_t>& colp2=sto2->colPointer();
  colp1.resize(nc+1);
  std::vector<number_t>::iterator it1=colp1.begin()+nc1;
  std::vector<number_t>::const_iterator it2;
  for (it2=colp2.begin(); it2!=colp2.end(); ++it2,++it1) *it1 = *it2 + nnz1;
  std::vector<number_t>& rowi1=sto1->rowIndex();
  const std::vector<number_t>& rowi2=sto2->rowIndex();
  rowi1.resize(nnz1+nnz2);
  it1=rowi1.begin()+nnz1;
  for (it2=rowi2.begin(); it2!=rowi2.end(); ++it2,++it1) *it1=*it2;
  sto1->nbOfColumns()=nc;

  //concatenate values
  if (vt==_real)  //real case
  {
    std::vector<real_t>& val1=mat1.rEntries_p->values();
    val1.resize(nnz1+nnz2+1);
    std::vector<real_t>& val2=mat2.rEntries_p->values();
    std::vector<real_t>::iterator it1=val1.begin()+nnz1+1;
    std::vector<real_t>::const_iterator it2;
    for (it2=val2.begin()+1; it2!=val2.end(); ++it2, ++it1) *it1=*it2;
  }
  else  //complex case
  {
    if (vt2==_real)
    {
      std::vector<complex_t>& val1 = mat1.cEntries_p->values();
      val1.resize(nnz1+nnz2+1);
      std::vector<real_t>& val2=mat2.rEntries_p->values();
      std::vector<complex_t>::iterator it1=val1.begin()+nnz1+1;
      std::vector<real_t>::const_iterator it2;
      for (it2=val2.begin()+1; it2!=val2.end(); ++it2, ++it1) *it1=*it2;
    }
    else
    {
      std::vector<complex_t>& val1=mat1.cEntries_p->values();
      val1.resize(nnz1+nnz2+1);
      std::vector<complex_t>& val2=mat2.cEntries_p->values();
      std::vector<complex_t>::iterator it1=val1.begin()+nnz1+1;
      std::vector<complex_t>::const_iterator it2;
      for (it2=val2.begin()+1; it2!=val2.end(); ++it2, ++it1) *it1=*it2;
    }
  }
  mat1.setNbOfCols(nc);

  //concatenate cdofsc
  cdofsc1.resize(nc1+nc2);
  std::vector<DofComponent>::iterator itc1=cdofsc1.begin()+nc1;
  std::vector<DofComponent>::const_iterator itc2=cdofsc2.begin();
  for (; itc2!=cdofsc2.end(); ++itc2, ++itc1) *itc1=*itc2;

  trace_p->pop();
}

// ---------------------------------------------------------------------------------------------------------------------------------
/*! merge constraints systems (if more than one conditions)
   it takes as input a vector of Constraints (each corresponding to an essential condition)
   and produces a map of Constraints indexed by unknown with two cases :

   Case of uncoupled unknowns (u1/v1, u2/v2 referred to same unknown u/v), u and v are not coupled by constraints
                 u1   u2   v1   v2
                --------------------                      u1   u2                       v1  v2
           c1   |cu1  0    0    0  |    | f1             ----------                    ---------
           c2   | 0  cu2   0    0  |  = | f2  ==>   Cu = |cu1  0  | = fu= | f1    Cv = |cv1 cv2| = fv= | f3
           c3   |     0   cv1  cv2 |    | f3             | 0  cu2 |       | f2         ---------
                --------------------                     ----------
    The merging process involving u1, u2 referring to same unknown u produces Constraints object where common dofs are merged
    One Constraints object for each unknown is created and returned

    Case of coupled unknowns (u1/v1, u2/v2 referred to same unknown u/v), u and v are coupled at least by one constraint
                u1    u2   v1   v2
                -------------------                       u   v
           c1   |cu1  0    0   0  |  = | f1             ---------
           c2   | 0  cu2  cv1 cv2 |    | f2     ==>     |Cu  Cv | = f
                -------------------                     ---------
    A global Constraints matrix is created and returned (indexed by 0)

    NOTE : when merging some Constraints in a new one the old ones are deleted by this function
*/
// ---------------------------------------------------------------------------------------------------------------------------------
std::map<const Unknown*, Constraints*> mergeConstraints(std::vector<Constraints*>& constraints)
{
  number_t nbc=constraints.size();
  if (nbc==0) error("is_null","constraints");
  trace_p->push("mergeConstraints");

  std::map<const Unknown*, Constraints*> mconstraints;
  std::vector<Constraints*>::iterator itc=constraints.begin();
  if (nbc==1)  //no merging
  {
    mconstraints[(*itc)->unknown()]= *itc;
    trace_p->pop();
    return mconstraints;
  }

  //collect constraints by unknown
  std::map<const Unknown*, std::list<Constraints*> > lsu;   //list of constraints involving same unknown
  std::map<const Unknown*, std::list<Constraints*> >::iterator itlsu;
  bool global = false;
  for (; itc!=constraints.end(); itc++)
  {
    std::set<const Unknown*> us=(*itc)->unknowns();
    std::set<const Unknown*>::iterator its=us.begin();
    for (; its!=us.end(); its++)
    {
      itlsu=lsu.find(*its);
      if (itlsu==lsu.end()) lsu[*its]=std::list<Constraints*>(1,*itc);
      else itlsu->second.push_back(*itc);
    }
    if (!global && (*itc)->coupledUnknowns()) global=true;
  }

  if (lsu.size()>1) global=true; // force global when more than one unknown involved

  std::set<Constraints*> undelete;   //set of constraints to not be deleted

  //merge constraints by unknown
  for (itlsu=lsu.begin(); itlsu!=lsu.end(); itlsu++)
  {
    const Unknown* u=itlsu->first;
    std::list<Constraints*>& ls=itlsu->second;
    std::list<Constraints*>::iterator itls;
    if (ls.size()==1 && !(*ls.begin())->coupledUnknowns())  //unique constraint on u and uncoupling constraint
    {
      mconstraints[u]=*ls.begin();
      undelete.insert(*ls.begin());   //no copy, do not delete after
    }
    else  //more than one constraints : merge them, if one coupling constraint : extract u part
    {
      //create new col cs storage and set up cdofs
      ValueType vtm=_real, vtr=_real;
      real_t vr; complex_t vc;
      bool newloc=true;
      bool newsym=true;
      EssentialConditions newconditions;
      std::map<DofComponent, std::set<number_t> > rowindex;
      std::map<DofComponent, std::set<number_t> >::iterator itmd, itmd2;
      std::vector<DofComponent>::iterator itd;
      std::vector<DofComponent> newcdofsc, newcdofsr;
      number_t nbr=0, nbc=0;
      for (itls=ls.begin(); itls!=ls.end(); itls++)
      {
        number_t c=1;
        for (itd=(*itls)->cdofsc_.begin(); itd!=(*itls)->cdofsc_.end(); itd++, c++) //loop on col cdofs
        {
          if (itd->u_p == u)
          {
            std::set<number_t> rls=(*itls)->matrixp()->storagep()-> getRows(c);  //row indices of col c
            std::set<number_t>::iterator its;
            if (nbr!=0) //add nbr to row indices (row translation)
            {
              std::set<number_t> srls;
              for (its=rls.begin(); its!=rls.end(); its++) srls.insert(*its + nbr);
              rls=srls;
            }
            if (rowindex.find(*itd)==rowindex.end()) rowindex[*itd]=rls;
            else rowindex[*itd].insert(rls.begin(), rls.end());
          }
        }
        newcdofsr.insert(newcdofsr.end(),(*itls)->cdofsr_.begin(),(*itls)->cdofsr_.end());  //update row cdofs (assuming no overlap in rows)
        nbr+=(*itls)->cdofsr_.size();
        //update data of merged constraint
        if ((*itls)->matrixp()->valueType_==_complex) vtm=_complex;  //update matrix value type
        if ((*itls)->rhsp()->valueType_==_complex) vtr=_complex;     //update rhs value type
        if (!(*itls)->local) newloc=false;                          //update local property
        if (!(*itls)->symmetric) newsym=false;                      //update symmetry property
        newconditions.insert(newconditions.end(),(*itls)->conditions_.begin(),(*itls)->conditions_.end());
      }
      nbr=newcdofsr.size();
      nbc=rowindex.size();
      newcdofsc.resize(nbc);
      itd= newcdofsc.begin();
      std::vector<std::vector<number_t> > vrowindex(nbc);
      std::vector<std::vector<number_t> >::iterator itv=vrowindex.begin();
      for (itmd=rowindex.begin(); itmd!=rowindex.end(); itmd++, itv++, itd++)
      {
        *itd = itmd->first;
        *itv = std::vector<number_t>(itmd->second.begin(),itmd->second.end());
      }
      MatrixStorage* ms=new ColCsStorage(nbr,nbc,vrowindex);
      //fill new constraint matrix and rhs vector
      MatrixEntry* newmatrix=new MatrixEntry(vtm, _scalar, ms);
      VectorEntry* newrhs= new VectorEntry(vtr,_scalar, nbr);
      nbr=0;
      for (itls=ls.begin(); itls!=ls.end(); itls++)
      {
        number_t c=1;
        const MatrixEntry* mec=(*itls)->matrixp();
        for (itd=(*itls)->cdofsc_.begin(); itd!=(*itls)->cdofsc_.end(); itd++, c++) //loop on col cdofs
        {
          if (itd->u_p == u)
          {
            itmd=rowindex.find(*itd);
            itmd2=rowindex.begin();
            number_t k=1;
            while (itmd2!=itmd) {k++; itmd2++;} //range of cdof in map
            std::set<number_t> rls=(*itls)->matrixp()->storagep()-> getRows(c);  //row indices of col c
            std::set<number_t>::iterator its;
            if (mec->valueType_==_real)
            {
              for (its=rls.begin(); its!=rls.end(); its++)
              {
                mec->getEntry(*its, c, vr);
                if (vtm==_real) newmatrix->setEntry(*its + nbr, k, vr);
                else newmatrix->setEntry(*its + nbr, k, complex_t(vr));
              }
            }
            else
            {
              for (its=rls.begin(); its!=rls.end(); its++)
              {
                mec->getEntry(*its, c, vc);
                newmatrix->setEntry(*its + nbr, k, vc);
              }
            }
          }
        }
        //fill new rhs
        const VectorEntry* vec=(*itls)->rhsp();
        if (vec->valueType_==_real)
        {
          for (number_t k=1; k<=vec->size(); k++)
          {
            vec->getEntry(k,vr);
            if (vtr==_real) newrhs->setEntry(k + nbr, vr);
            else newrhs->setEntry(k + nbr, complex_t(vr));
          }
        }
        else
        {
          for (number_t k=1; k<=vec->size(); k++)
          {
            vec->getEntry(k,vc);
            newrhs->setEntry(k + nbr, vc);
          }
        }
        nbr+=(*itls)->cdofsr_.size();
      }
      //create new constraint
      Constraints* newc= new Constraints(newmatrix,newrhs);
      newc->cdofsc_=newcdofsc;
      newc->cdofsr_=newcdofsr;
      newc->conditions_=newconditions;
      newc->local = newloc;
      newc->symmetric=newsym;
      mconstraints[u]=newc;
    }
  }

  //global merge if some constraints couple different unknowns
  //merge all previous constraints in one constraint
  if (global)
  {
    std::map<const Unknown*, Constraints*>::iterator itmc;
    std::map<DofComponent,number_t> rowdofs;
    std::map<DofComponent,number_t>::iterator itmn;
    std::vector<DofComponent> newcdofsc, newcdofsr;
    std::vector<DofComponent>::iterator itd;
    ValueType vtm=_real, vtr=_real;
    number_t k=1;
    for (itmc=mconstraints.begin(); itmc!=mconstraints.end(); itmc++)
    {
      Constraints* cons=itmc->second;
      for (itd=cons->cdofsr_.begin(); itd!=cons->cdofsr_.end(); itd++)
        if (rowdofs.find(*itd)==rowdofs.end())
          {rowdofs.insert(std::make_pair(*itd,k)); k++;}
      newcdofsc.insert(newcdofsc.end(),cons->cdofsc_.begin(),cons->cdofsc_.end());
    }
    newcdofsr.resize(rowdofs.size());
    itd=newcdofsr.begin();
    for (itmn=rowdofs.begin(); itmn!=rowdofs.end(); itmn++, itd++) *itd = itmn->first;
    number_t nbc=newcdofsc.size(), nbr= newcdofsr.size();
    std::vector<std::vector<number_t> > vrowindex(nbc);
    std::vector<std::vector<number_t> >::iterator itv=vrowindex.begin();
    std::set<number_t>::iterator its;
    std::vector<number_t>::iterator itvr;
    for (itmc=mconstraints.begin(); itmc!=mconstraints.end(); itmc++)
    {
      Constraints* cons=itmc->second;
      number_t c=1;
      for (itd=cons->cdofsc_.begin(); itd!=cons->cdofsc_.end(); itd++, itv++, c++)
      {
        std::set<number_t> rls=cons->matrixp()->storagep()-> getRows(c);  //row indices of col c
        std::vector<number_t> vrls(rls.size());
        itvr=vrls.begin();
        for (its=rls.begin(); its!=rls.end(); its++, itvr++) *itvr = rowdofs[cons->cdofsr_[*its-1]];
        *itv = vrls;
        if (cons->matrixp()->valueType_ == _complex) vtm = _complex;
        if (cons->rhsp()->valueType_ == _complex) vtr = _complex;
      }
    }
    MatrixStorage* ms=new ColCsStorage(nbr,nbc,vrowindex);   //new storage
    MatrixEntry* newmatrix=new MatrixEntry(vtm, _scalar, ms);
    VectorEntry* newrhs= new VectorEntry(vtr,_scalar, nbr);
    k=1; real_t vr; complex_t vc;
    for (itmc=mconstraints.begin(); itmc!=mconstraints.end(); itmc++)
    {
      Constraints* cons=itmc->second;
      //update matrix
      number_t c=1;
      for (itd=cons->cdofsc_.begin(); itd!=cons->cdofsc_.end(); itd++, itv++, k++, c++)
      {
        std::set<number_t> rls=cons->matrixp()->storagep()-> getRows(c);
        if (cons->matrixp()->valueType_ ==_real)
        {
          for (its=rls.begin(); its!=rls.end(); its++)
          {
            number_t r = rowdofs[cons->cdofsr_[*its-1]];
            cons->matrixp()->getEntry(*its, c, vr);
            if (vtm==_real) newmatrix->setEntry(r,k,vr);
            else newmatrix->setEntry(r,k,complex_t(vr));
          }
        }
        else
        {
          for (its=rls.begin(); its!=rls.end(); its++, itvr++)
          {
            number_t r = rowdofs[cons->cdofsr_[*its-1]];
            cons->matrixp()->getEntry(*its, c, vc);
            newmatrix->setEntry(r,k,vc);
          }
        }
      }
      //update rhs
      const VectorEntry* vec=cons->rhsp();
      number_t k=1;
      for (itd=cons->cdofsr_.begin(); itd!=cons->cdofsr_.end(); itd++, k++)
      {
        number_t r = rowdofs[*itd];
        if (vec->valueType_==_real)
        {
          vec->getEntry(k,vr);
          if (vtr==_real) newrhs->setEntry(r, vr);
          else newrhs->setEntry(r, complex_t(vr));
        }
        else
        {
          vec->getEntry(k,vc);
          newrhs->setEntry(r, vc);
        }
      }
    }
    //delete mconstraints and reallocate
    for (itmc=mconstraints.begin(); itmc!=mconstraints.end(); itmc++)
    {
      if (undelete.find(itmc->second) == undelete.end())
      {
        undelete.insert(itmc->second);
        delete itmc->second;
      }
    }
    mconstraints.clear();
    Constraints* newc= new Constraints(newmatrix,newrhs);
    newc->cdofsc_=newcdofsc;
    newc->cdofsr_=newcdofsr;
    bool newloc=true;
    bool newsym=true;
    EssentialConditions newconditions;
    for (itc=constraints.begin(); itc!=constraints.end(); itc++)
    {
      newconditions.insert(newconditions.end(),(*itc)->conditions_.begin(),(*itc)->conditions_.end());
      if (!(*itc)->local) newloc=false;
      if (!(*itc)->symmetric) newsym=false;
    }
    newc->conditions_= newconditions;
    newc->local = newloc;
    newc->symmetric = newsym;
    mconstraints[0]=newc;
  }

  //delete old constraints
  for (itc=constraints.begin(); itc!=constraints.end(); itc++)
    if (undelete.find(*itc) == undelete.end()) delete *itc;

  //end of merging
  trace_p->pop();
  return mconstraints;
}


// ---------------------------------------------------------------------------------------------------------------------------------
/*! reduce constraints system to an upper triangular system using stable QR factorization
    transform the constraints into a new constraints where matrix is an upper triangular matrix

                           u1           u2         ...      un
                     --------------------------------------------
                     | ... ... ... | ... ... |           |  ... |          | . |
    original     C = | ...  C1 ... | ... ... |    ...    |   Cn |       D= | . |
    constraints      | ... ... ... | ... ... |           |  ... |          | . |
                     --------------------------------------------

                           rearrangement of u1 , u2, ..., un
                     --------------------------------------------
                     | 1  x ... | ... ... |           |  ... ...|          | . |
    reduce      rC = | 0  1  x  | ... ... |    ...    |  ... ...|      rD= | . |
    constraints      | 0  0 ... | ... ... |           |  ... ...|          | . |
                     --------------------------------------------

  when algorithm 'fails' it means that some constraints are redundant or conflicting, so we eliminate on the fly some line constraints

  In a second step, the upper triangular system is inverted to produce residual rectangular matrix and modified right hand side

                                    part of  u1 , u2, ..., un
                                --------------------------------
                                |  x   x  |           |  x   x |          | . |
    residual constraints   F =  |  x  ... |    ...    |  x  ...|      b = | . |
                                | ... ... |           | ... ...|          | . |
                                --------------------------------
    The first p column indices  (of squared upper triangular matrix) corresponds to the eliminated dofs (say U_e)
    while the other column indices corresponds to the residual dofs (say U_r). Finally the constraints system reads

                            ----------------------------------------
                            | U_e + F*U_r = b  <=>  U_e= b - F*U_r |
                            ----------------------------------------

    which is the useful expression to deal with essential conditions in bilinear form computation
    \note F may be a null matrix (classic Dirichlet condition for instance, condition reads U_e = b as usual)

  \param aszero  : a positive value near 0 used to round to 0 very small values, say |v|<aszero, default value is 10*theEpsilon

  this process may be permutes column cdofs and produces
  elcdofs_ : the map of eliminated cdofs, cdof -> k (rank in cdofsc_)
  recdofs_ : the map of reduced cdofs : cdof -> k (rank in cdofsc_)

*/
// ---------------------------------------------------------------------------------------------------------------------------------

void Constraints::reduceConstraints(real_t aszero)
{
  trace_p->push("Constraints::reduceConstraints");

  if (isId)  //special case of Id : trivial reduction
  {
    //set elcdofs and recdofs
    number_t n=matrix_p->nbOfRows();
    std::vector<DofComponent>::iterator itc=cdofsc_.begin();
    for (number_t k=1; k<=n; k++, itc++) elcdofs_[*itc]=k;
    recdofs_.clear();
    // matrix_p->clear();
    // delete matrix_p;
    // matrix_p=0;
    if (rhs_p!=0 && norminfty(*rhs_p)< aszero)  {delete rhs_p; rhs_p=0;}  //delete near 0 right hand side
    trace_p->pop();
    return;
  }

  //General case, use QR reduction
  MatrixEntry matR, matQ;   //for QR factorisation
  bool withPerm = true;     //column permutation allowed in QR factorisation
  std::vector<number_t>* numcol=0;
  number_t stop=0;
  QR(*matrix_p, matR, false, matQ, rhs_p, withPerm, numcol, stop);
  number_t nbr=matrix_p->nbOfRows(), nbc=matrix_p->nbOfCols();
  // thePrintStream << "matR:" << eol << matR << eol;
  // thePrintStream << "stop=" << stop << " numcol=" << *numcol << eol;

  if (stop < nbr)
  {
    VectorEntry crhs(*rhs_p);
    crhs.deleteRows(1,stop);
    number_t nbz=crhs.nbzero(aszero), nbnz=nbr-stop-nbz;
    std::stringstream ss;
    ss << "Constraints::reduceConstraints() : in essential conditions \n" << conditions_;
    if (nbnz>0) ss << " conflicting constraint row(s) detected, up to " << (nbr-stop) << " rows" << eol;
    else  ss << nbz << " redundant constraint row(s) detected and eliminated" << eol;
    warning("free_warning",ss.str());
    matR.deleteRows(stop+1, nbr);  //delete redundant rows
    rhs_p->deleteRows(stop+1, nbr);
    cdofsr_.resize(stop);          //update cdofsr_
  }
  if (withPerm) //apply permutation to cdofsc_
  {
    std::vector<DofComponent> newdof(cdofsc_.size());
    std::vector<DofComponent>::iterator itc = newdof.begin();
    std::vector<number_t>::iterator itn=numcol->begin();
    for (; itn!=numcol->end(); itn++, itc++) *itc = cdofsc_[*itn];
    cdofsc_=newdof;
  }

  //set elcdofs and recdofs
  nbr=matR.nbOfRows();
  std::vector<DofComponent>::iterator itc=cdofsc_.begin();
  number_t k=1;
  for (; k<=nbr; k++, itc++) elcdofs_[*itc]=k;
  for (; k<=nbc; k++, itc++) recdofs_[*itc]=k;

  // reduced upper triangular system to "diagonal" system
  matR.roundToZero(aszero);     //to eliminate spurious rounding effects
  MatrixEntry matU(matR,true);  //copy matR, forcing storage copy to prevent shared storage
  nbc=matR.nbOfCols();
  if (nbc>nbr)
  {
    matU.deleteCols(nbr+1,nbc);
    matR.deleteCols(1, nbr);
    if (!matU.isId(theTolerance)) QRSolve(matU,&matR,rhs_p);  // solve upper triangular system

    if(matR.norminfty()<theEpsilon) //clear reduced matrix if null matrix
    {
       matR.clear();
       recdofs_.clear(); //no reduced dofs
    }
    else
      *matrix_p=matR;  //ERIC : check the storage management, delete, copy, shared storage ...
  }
  else
  {
    matR.clear();
    // delete matrix_p;
    // matrix_p=0;
    if (!matU.isId(theTolerance)) QRSolve(matU,0,rhs_p);     // solve upper triangular system
  }

  //finalization
  if (rhs_p!=0 && norminfty(*rhs_p)< aszero)  //delete near 0 right hand side
  {
    delete rhs_p;
    rhs_p=0;
  }
  if (numcol!=0) delete numcol;
  reduced=true;
  trace_p->pop();
}

// ---------------------------------------------------------------------------------------------------------------------------------
/*! extract eliminated and reduced cdofs from list of cdofs
   cdofs : list of cdofs
   melcdofs : list of elimited cdofs with their ranks in cdofs list
   mrecdofs : list of reduced cdofs with their ranks in cdofs list
   useDual  : if true use either unknown or dual unknown to get cdofs
*/
// ---------------------------------------------------------------------------------------------------------------------------------
void Constraints::extractCdofs(const std::vector<DofComponent>& cdofs, std::map<DofComponent, number_t>& melcdofs,
                               std::map<DofComponent, number_t>& mrecdofs, bool useDual) const
{
  melcdofs.clear(); mrecdofs.clear();
  std::map<DofComponent, number_t>::const_iterator itmd;
  std::vector<DofComponent>::const_iterator itd;
  std::vector<DofComponent>::iterator itdd;
  std::vector<DofComponent> dual_cdofs;
  number_t k=1;
  if (useDual)
  {
    dual_cdofs.resize(cdofs.size());
    std::vector<DofComponent>::iterator itdd = dual_cdofs.begin();
    for (itd=cdofs.begin(); itd!=cdofs.end(); itd++, itdd++) *itdd = itd->dual();
  }
  itdd=dual_cdofs.begin();
  for (itd=cdofs.begin(); itd!=cdofs.end(); itd++,k++) //travel all row cdof of matrix mat
  {
    itmd=elcdofs_.find(*itd);
    if (itmd!=elcdofs_.end()) melcdofs[*itd]=k;     // eliminated cdof found
    itmd=recdofs_.find(*itd);
    if (itmd!=recdofs_.end()) mrecdofs[*itd]=k;     // reduced cdof found
    if (useDual)
    {
      itmd=elcdofs_.find(*itdd);
      if (itmd!=elcdofs_.end()) melcdofs[*itd]=k;     // eliminated cdof found
      itmd=recdofs_.find(*itdd);
      if (itmd!=recdofs_.end()) mrecdofs[*itd]=k;     // reduced cdof found
      itdd++;
    }
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------
// print utilities
// ---------------------------------------------------------------------------------------------------------------------------------
void Constraints::print(std::ostream& os) const
{
  os << " Constraints system related to essential condition(s)" << eol << conditions_;
  if (matrix_p==0)
  {
    os << " no matrix representation !" << eol;
    return;
  }
  if (coupledUnknowns()) os << " constraints system couples different unknowns, global system ";
  //else os << " constraints system on unknown " << unknown()->name()  <<  " ";
  if (unknown()!=0) os << " constraints system on unknown " << unknown()->name()  <<  " ";
  else
  {
    std::set<const Unknown*> us=conditions_.unknowns();
    std::set<const Unknown*>::iterator itu = us.begin();
    if (isTestMode)
    {
      std::set<string_t> sus;
      for (; itu!=us.end(); ++itu) sus.insert((*itu)->name());
      std::set<string_t>::iterator itsu = sus.begin();
      os << " constraints system on unknowns (" << (*itsu);
      ++itsu;
      for (; itsu!=sus.end(); ++itsu) os << ", " << (*itsu);
    }
    else
    {
      os << " constraints system on unknowns (" << (*itu)->name(); ++itu;
      for (; itu!=us.end(); ++itu) os << ", " << (*itu)->name();
    }
    os << ") ";
  }
  if (rhs_p==0) os << " (homogeneous)";
  else os << " (non homogeneous)";
  if (local) os << " , local constraints";
  else os << " , non local constraints";
  if (symmetric) os << " , keeping symmetry";
  else os << " , not keeping symmetry";
  os << eol;
  os << "  unknown dofs involved (" << cdofsc_.size() << ") : " <<  cdofsc_ << eol;
  os << "  test function dofs involved (" << cdofsr_.size() << ") : " <<  cdofsr_ << eol;
  os << "  constraints matrix : " << *matrix_p;
  if (rhs_p!=0) os << "  right hand side : " << *rhs_p << eol;
  if (reduced)
  {
    os << "  eliminated dofs (" << elcdofs_.size() << ") : " <<  elcdofs_ << eol;
    os << "  reduced dofs (" << recdofs_.size() << ") : " <<  recdofs_ << eol;
  }
  return;
}

std::ostream& operator<<(std::ostream& os, const Constraints& cd)
{
  cd.print(os);
  return os;
}

//utility to create a unique unknown name associated to an essential condition
string_t unknownEcName(const EssentialCondition& ec)
{
  string_t nar=ec.name()+" (c";
  number_t k=0;
  string_t na=nar+tostring(k)+")";
  while (findUnknown(na)!=0) { k++; na=nar+tostring(k)+")";}
  return na;
}

/*!
  perform pseudo reduction of reduced essential conditions in a matrix. Essential conditions have the following form  :
                          U_E + F*U_R = H   for column unknown U
                          V_E + G*V_R = 0   for row test function V (generally related to unknown U)
  where E are the eliminated unknown/test function indices and R are the reduced unknown/test function indices
  other indices S correspond to dofs not connected to the essential conditions dofs
  The pseudo reduction of matrix consists in
    - modifying column A.j for j in R by adding -Fkj*A.k for k in E and replacing column A.k for k in E by a 0 column
    - modifying row Ai. for i in R by adding -Gki*Ak. for k in E and replacing row Ak. for k in E by a 0 row
    - if eliminated v unknowns are dual of eliminated u unknowns (G=F), the VE rows are replaced by the equation (or a part of)
                                      a*U_E + a*F*U_R = a*H   where a is a given scalar
  to delay the right hand side modification, the (A.k) columns (k in E) are stored in a new structure

  At the end of the process, the eliminated system looks like (C the correction matrix mxE get from reduction of the constaints)

                     U_E       U_R     U_S        U        RHS
                 ----------------------------   -------   -------
                 |          |        |      |   |     |   |     |
            V_E  | (a+b)*Id |  a*F   |   0  |   | U_E |   | a*H |    => (a+b)*U_E + a*F*U_R = a*H (mimics the constraints)
                 |          |        |      |   |     |   |     |
                 ----------------------------   -------   -------
                 |          |        |      |   |     |   |     |
            V_R  |     0    |  A_RR  | A_RS | * | U_R | = | B_R |    = B0_R - C_RE*H -Gt*B0_E  (rhs correction on reduced dof)
                 |          |        |      |   |     |   |     |
                 ----------------------------   -------   ------
                 |          |       |       |   |     |   |     |
            V_S  |     0    |  A_SR | A_SS  |   | U_S |   | B_S |    = B0_S - C_SE*H           (rhs correction on free dofs)
                 |          |       |       |   |     |   |     |
                 ----------------------------   -------   -------

  If F=G=0 (Dirichlet condition) , R={} and the system reads

                     U_E       U_S         U        RHS
                 --------------------   -------   -------
                 |          |       |   |     |   |     |
            V_E  | (a+b)*Id |   0   |   | U_E |   | a*H |    => (a+b)*U_E = a*H  (constraints)
                 |          |       |   |     |   |     |
                 -------------------  * ------- = -------
                 |          |       |   |     |   |     |
            V_S  |     0    |  A_SS |   | U_S |   | B_S |    = B0_S - C_SE*H    (rhs correction)
                 |          |       |   |     |   |     |
                 --------------------   -------   -------
                 for homogeneous dirichlet condition H=0 and B_S = B0_S!

  In some cases (F=G=0 or non dof coupling condition...) the storage is not modified but in other cases (transmission condition for instance)
  the storage is modified

  choosing a=0,b!=0 allows to deal with eigen value problems by shifing the spectra corresponding to eliminated dof

    mat   : matrix to be eliminated
    cdofr : row component dofs of matrix
    cdofc : col component dofs of matrix
    rhsmat: right hand side matrix pointer

*/
void Constraints::pseudoColReduction(MatrixEntry* mat, std::vector<DofComponent>& cdofr,
                                     std::vector<DofComponent>& cdofc, MatrixEntry*& rhsmat)
{
  trace_p->push("Constraints::pseudoColReduction()");
  if (matrix_p==0) {trace_p->pop(); return;}       //no constraints to apply

  //locate eliminated cdofs and reduced cdofs in matrix columns
  std::map<DofComponent, number_t> melcdofs, mrecdofs;
  extractCdofs(cdofc, melcdofs, mrecdofs, false);
  number_t nbe=elcdofs_.size();
  std::map<DofComponent, number_t>::iterator itn, itm, itmd;
  std::vector<DofComponent>::iterator itd;
  complex_t z0(0.,0.);

  if (melcdofs.size()>0) // create correction matrix (rshmat)
  {
    //construct storage and constraint correction matrix
    MatrixStorage* ms=0;
    if (mat->storageType()!=_dense)
    {
      std::vector< std::vector<number_t> > rowindex(melcdofs.size());
      std::vector< std::vector<number_t> >::iterator itri=rowindex.begin();
      for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++, itri++)
      {
        std::set<number_t> snum=mat->storagep()->getRows(itn->second);
        *itri = std::vector<number_t>(snum.begin(), snum.end());
      }
      ms= new ColCsStorage(mat->nbOfRows(),melcdofs.size(),rowindex);
    }
    else //dense storage
      ms= new ColDenseStorage(mat->nbOfRows(),melcdofs.size());

    rhsmat = new MatrixEntry(mat->valueType_, _scalar, ms);
    number_t j=1;
    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++, itd++, j++)
    {
      number_t c=itn->second;
      std::vector<std::pair<number_t, number_t> > adrs = mat->storagep()->getCol(mat->symmetry(),c);
      std::vector<number_t> mat_adr(adrs.size());
      std::vector<std::pair<number_t, number_t> >::iterator ita=adrs.begin();
      std::vector<number_t>::iterator itma=mat_adr.begin();
      for (; ita!=adrs.end(); ita++, itma++) *itma = ita->second;   //addresses in mat of col c
      adrs=rhsmat->storagep()->getCol(_noSymmetry,j);
      std::vector<number_t> rhsmat_adr(adrs.size());
      itma=rhsmat_adr.begin();
      for (ita=adrs.begin(); ita!=adrs.end(); ita++, itma++) *itma = ita->second;  //addresses in rhsmat of col c
      rhsmat->copyVal(*mat, mat_adr, rhsmat_adr);
    }
  }

  if (mrecdofs.size()>0) // update reduced columns
  {
    //change mat type if  real matrix and complex constraints matrix
    ValueType vtm=mat->valueType_, vtcm=matrix_p->valueType_;
    if (vtm==_real && vtcm==_complex)
    {
      mat->toComplex();
      mat->valueType_ = _complex;
      vtm=_complex;
    }
    // dual or not dual unknown
    bool meldual = false, recdual=false;
    if(melcdofs.size()>0) meldual = (elcdofs_.find(melcdofs.begin()->first) == elcdofs_.end());
    if(mrecdofs.size()>0) recdual = (recdofs_.find(mrecdofs.begin()->first) == recdofs_.end());

    real_t cr=0.; complex_t cc=z0;
    //combine columns
//    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
//    {
//      std::map<DofComponent, number_t>::const_iterator itmd;
//      if(meldual) itmd = elcdofs_.find(itn->first.dual()); else itmd =elcdofs_.find(itn->first);
//      number_t i=itmd->second, c=itn->second;
//      for (itm=mrecdofs.begin(); itm!=mrecdofs.end(); itm++)
//      {
//        if(recdual) itmd = recdofs_.find(itm->first.dual()); else itmd =recdofs_.find(itm->first);
//        number_t j=itmd->second - nbe;
//        if (vtcm==_real)
//        {
//          matrix_p->getEntry(i,j,cr);
//          if (cr!=0)
//          {
//            if (vtm == _real) mat->addColToCol(c,itm->second, -cr);
//            else mat->addColToCol(c,itm->second,complex_t(-cr));
//          }
//        }
//        else
//        {
//          matrix_p->getEntry(i,j,cc);
//          if (cc!=z0) mat->addColToCol(c,itm->second,-cc);
//        }
//      }
//    }
    //reverse loop to avoid data race in omp
    std::vector<number_t> melind(melcdofs.size());  // move map into vectors to avoid find in internal loop
    std::vector<number_t> melpos(melcdofs.size());
    std::vector<number_t>::iterator iti=melind.begin(), itp=melpos.begin();
    for (std::map<DofComponent, number_t>::iterator itn=melcdofs.begin(); itn!=melcdofs.end(); ++itn, ++itp, ++iti)
    {
      *itp = itn->second;
      if (meldual) *iti = elcdofs_.find(itn->first.dual())->second;
      else *iti=elcdofs_.find(itn->first)->second;
    }
    #ifdef XLIFEPP_WITH_OMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(number_t n=0;n<mrecdofs.size();n++)
    {
      std::map<DofComponent, number_t>::iterator itm=mrecdofs.begin();
      std::advance(itm,n);
      number_t c=itm->second, j;
      if (recdual) j = recdofs_.find(itm->first.dual())->second - nbe;
      else j = recdofs_.find(itm->first)->second - nbe;
      std::vector<number_t>::iterator itp = melpos.begin();
      for (std::vector<number_t>::iterator iti=melind.begin(); iti!=melind.end();++iti,++itp)
      {
        real_t cr; complex_t cc;
        if (vtcm==_real)
        {
          matrix_p->getEntry(*iti,j,cr);
          if (cr!=0)
          {
            if (vtm == _real) mat->addColToCol(*itp,c, -cr);
            else mat->addColToCol(*itp,c,complex_t(-cr));
          }
        }
        else
        {
          matrix_p->getEntry(*iti,j,cc);
          if (cc!=z0) mat->addColToCol(*itp,c,-cc);
        }
      }
    }
  }

  if (melcdofs.size()>0) // eliminate columns
  {
    //pseudo col reduction (set eliminated columns to 0)
    //for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++) mat->setColToZero(itn->second,itn->second);
    #ifdef XLIFEPP_WITH_OMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(number_t n=0;n<melcdofs.size();n++)
    {
      std::map<DofComponent, number_t>::iterator itn=melcdofs.begin();
      std::advance(itn,n);
      number_t c=itn->second;
      mat->setColToZero(c,c);
    }
  }
  trace_p->pop();
}



// pseudo reduction of rows
//   mat : matrix to be reduced
//   cdofr : row dofs
//   cdofc : column dofs
//   alpha, beta : reduction coefficients
//   subConstraint : if true, substitute the eliminated x reduced block by the constraints matrix
void Constraints::pseudoRowReduction(MatrixEntry* mat, std::vector<DofComponent>& cdofr, std::vector<DofComponent>& cdofc,
                                     const complex_t& alpha, const complex_t& beta, bool subConstraint)
{
  trace_p->push("Constraints::pseudoRowReduction()");
  if (matrix_p==0) {trace_p->pop(); return;}       //no constraints to apply

  //locate eliminated cdofs and reduced cdofs
  std::map<DofComponent, number_t> melcdofs, mrecdofs;
  extractCdofs(cdofr, melcdofs, mrecdofs, true);
  number_t nbe=elcdofs_.size();
  ValueType vtm=mat->valueType_, vtcm=matrix_p->valueType_;
  complex_t z0(0.,0.);
  real_t cr=0.; complex_t cc=z0;

  std::map<DofComponent, number_t>::iterator itn, itm, itmd;
  std::vector<DofComponent>::iterator itd;
  //elapsedTime("enter row combination",theCout);
  if (mrecdofs.size()>0) // update reduced rows
  {
    //change mat type if required (real matrix and complex constraints matrix)
    if (vtm==_real && vtcm==_complex)
    {
      mat->toComplex();
      mat->valueType_ = _complex;
      vtm=_complex;
    }
    // dual or not dual unknown
    bool meldual = false, recdual=false;
    if(melcdofs.size()>0) meldual = (elcdofs_.find(melcdofs.begin()->first) == elcdofs_.end());
    if(mrecdofs.size()>0) recdual = (recdofs_.find(mrecdofs.begin()->first) == recdofs_.end());
    //combine rows
//    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
//    {
//        std::map<DofComponent, number_t>::const_iterator itmd;
//        if (meldual) itmd=elcdofs_.find(itn->first.dual()); else itmd=elcdofs_.find(itn->first);
//        number_t i=itmd->second;
////      for (itm=mrecdofs.begin(); itm!=mrecdofs.end(); itm++)
////      {
//      for(number_t n=0;n<mrecdofs.size();n++)
//      {
//        std::map<DofComponent, number_t>::iterator itm=mrecdofs.begin();
//        std::advance(itm,n);
//        std::map<DofComponent, number_t>::const_iterator itmdr;
//        if (recdual) itmdr=recdofs_.find(itm->first.dual()); else itmdr=recdofs_.find(itm->first);
//        number_t j=itmdr->second - nbe;
//        real_t cr; complex_t cc;
//        if (vtcm==_real)
//        {
//          matrix_p->getEntry(i,j,cr);
//          if (cr!=0)
//          {
//            if (vtm == _real) mat->addRowToRow(itn->second,itm->second, -cr);
//            else mat->addRowToRow(itn->second,itm->second,complex_t(-cr));
//          }
//        }
//        else
//        {
//          matrix_p->getEntry(i,j,cc);
//          if (cc!=z0) mat->addRowToRow(itn->second,itm->second,conj(-cc));  //warning : conjugate constraints coefficient
//        }
//      }
//    }
//  }
      //reverse loop to avoid data race
      std::vector<number_t> melind(melcdofs.size());  // move map into vectors to avoid find in internal loop
      std::vector<number_t> melpos(melcdofs.size());
      std::vector<number_t>::iterator iti=melind.begin(), itp=melpos.begin();
      for (std::map<DofComponent, number_t>::iterator itn=melcdofs.begin(); itn!=melcdofs.end(); ++itn, ++itp, ++iti)
      {
        *itp = itn->second;
        if (meldual) *iti = elcdofs_.find(itn->first.dual())->second;
        else *iti=elcdofs_.find(itn->first)->second;
      }
    //elapsedTime("row combination preparation",theCout);
    #ifdef XLIFEPP_WITH_OMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(number_t n=0;n<mrecdofs.size();n++)
    {
      std::map<DofComponent, number_t>::iterator itm=mrecdofs.begin();
      std::advance(itm,n);
      number_t r=itm->second, j;
      if (recdual) j = recdofs_.find(itm->first.dual())->second - nbe;
      else j = recdofs_.find(itm->first)->second - nbe;
      std::vector<number_t>::iterator itp = melpos.begin();
      for (std::vector<number_t>::iterator iti=melind.begin(); iti!=melind.end();++iti,++itp)
      {
        real_t cr; complex_t cc;
        if (vtcm==_real)
        {
          matrix_p->getEntry(*iti,j,cr);
          if (cr!=0)
          {
            if (vtm == _real) mat->addRowToRow(*itp,r, -cr);
            else mat->addRowToRow(*itp,r,complex_t(-cr));
          }
        }
        else
        {
          matrix_p->getEntry(*iti,j,cc);
          if (cc!=z0) mat->addRowToRow(*itp,r,conj(-cc));  //warning : conjugate constraints coefficient
        }
      }
    }
  }
  //elapsedTime("row combination",theCout);

  if (melcdofs.size()>0) // "eliminate" rows
  {
//    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
//        mat->setRowToZero(itn->second,itn->second);
    #ifdef XLIFEPP_WITH_OMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(number_t n=0;n<melcdofs.size();n++)
    {
      std::map<DofComponent, number_t>::iterator itn=melcdofs.begin();
      std::advance(itn,n);
      number_t r=itn->second;
      mat->setRowToZero(r,r);
    }
    //elapsedTime("set row to zero",theCout);
    //set diagonal coefficient of mat to alpha+beta if row cdofs and col cdofs are dual or same
    number_t k=1, kt;
    std::map<DofComponent, number_t> mapcdofc;  //map of column dofs
    for (itd=cdofc.begin(); itd!=cdofc.end(); itd++, k++) mapcdofc[*itd]=k;
    ValueType vtm=mat->valueType_;
    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
    {
      k=itn->second;
      itm=mapcdofc.find(itn->first);
      if (itm==mapcdofc.end()) itm=mapcdofc.find(itn->first.dual());
      if (itm!=mapcdofc.end())
      {
        kt = itm->second;
        if (vtm==_real) mat->setEntry(k,kt, alpha.real()+beta.real());
        else  mat->setEntry(k,kt,alpha+beta);
      }
    }
  }
  //elapsedTime("eliminate rows",theCout);

  if (subConstraint && mrecdofs.size()>0) // replace "eliminate" row assuming v_constraints is same as u_constraints
  {
    //replace block eliminated x reduced by constraints matrix
    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
    {
      std::map<DofComponent, number_t>::const_iterator itmd=elcdofs_.find(itn->first);
      if (itmd==elcdofs_.end()) itmd=elcdofs_.find(itn->first.dual());  //try with dual
      number_t i=itmd->second;
      for (itm=mrecdofs.begin(); itm!=mrecdofs.end(); itm++)
      {
        itmd=recdofs_.find(itm->first);
        if (itmd==recdofs_.end()) itmd=recdofs_.find(itm->first.dual());  //try with dual

        number_t j=itmd->second - nbe; //block eliminated x reduced
        if (vtcm==_real)
        {
          matrix_p->getEntry(i,j,cr);
          if (cr!=0.)  //all rows have been already reset to 0 (except diagonal to alpha)
          {
            if (vtm==_real) mat->setEntry(itn->second, itm->second,alpha.real()*cr,false);
            else mat->setEntry(itn->second, itm->second,alpha*complex_t(cr),false);
          }
        }
        else
        {
          matrix_p->getEntry(i,j,cc);
          if (cc!=z0)
          {
            mat->setEntry(itn->second, itm->second,alpha*cc,false);
          }
        }
      }
    }
  }
  //elapsedTime("insert constraints",theCout);
  trace_p->pop();
}


/*! penalization reduction
    add alpha(U_E + F*U_R) to the block (Ve,Ue+Ur) of the matrix mat, assuming same constraints on u and v
    works only for Dirichlet condition (no reduced unknown component)
*/
void Constraints::penalizationReduction(MatrixEntry* mat, std::vector<DofComponent>& cdofr, std::vector<DofComponent>& cdofc,
                                        const complex_t& alpha)
{
  trace_p->push("Constraints::penalizationReduction()");
  if (matrix_p==0) {trace_p->pop(); return;}       //no constraints to apply

//locate eliminated cdofs and reduced cdofs
  std::map<DofComponent, number_t> melcdofs, mrecdofs;
  extractCdofs(cdofr, melcdofs, mrecdofs, true);
  number_t nbe=elcdofs_.size();
  ValueType vtm=mat->valueType_, vtcm=matrix_p->valueType_;
  complex_t z0(0.,0.);
  real_t cr=0., cmr=0.; complex_t cc=z0, cmc=z0;

  if (mrecdofs.size()>0) warning("free_warning","penalization of non Dirichlet condition is hazardous, use pseudo-reduction");

  std::map<DofComponent, number_t>::iterator itn, itm, itmd;
  std::vector<DofComponent>::iterator itd;

  if (melcdofs.size()>0) //add alpha to the  diagonal coefficient of mat if row cdofs and col cdofs are dual or same
  {
    number_t k=1, kt;
    std::map<DofComponent, number_t> mapcdofc;  //map of column dofs
    for (itd=cdofc.begin(); itd!=cdofc.end(); itd++, k++) mapcdofc[*itd]=k;
    ValueType vtm=mat->valueType_;
    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
    {
      k=itn->second;
      itm=mapcdofc.find(itn->first);
      if (itm==mapcdofc.end()) itm=mapcdofc.find(itn->first.dual());
      if (itm!=mapcdofc.end())
      {
        kt = itm->second;
        if (vtm==_real)
        {
          mat->getEntry(k,kt,cr);
          mat->setEntry(k,kt, cr+alpha.real());
        }
        else
        {
          mat->getEntry(k,kt,cc);
          mat->setEntry(k,kt,cc+alpha);
        }
      }
    }
  }

  if (mrecdofs.size()>0) // assuming v_constraints is same as u_constraints
  {
    //add constraints to eliminated x reduced block
    for (itn=melcdofs.begin(); itn!=melcdofs.end(); itn++)
    {
      std::map<DofComponent, number_t>::const_iterator itmd=elcdofs_.find(itn->first);
      if (itmd==elcdofs_.end()) itmd=elcdofs_.find(itn->first.dual());  //try with dual
      number_t i=itmd->second;
      for (itm=mrecdofs.begin(); itm!=mrecdofs.end(); itm++)
      {
        itmd=recdofs_.find(itm->first);
        if (itmd==recdofs_.end()) itmd=recdofs_.find(itm->first.dual());  //try with dual

        number_t j=itmd->second - nbe; //block eliminated x reduced
        if (vtcm==_real)
        {
          matrix_p->getEntry(i,j,cmr);
          if (cmr!=0.)
          {
            if (vtm==_real)
            {
              mat->getEntry(itn->second, itm->second,cr);
              mat->setEntry(itn->second, itm->second,cr+alpha.real()*cmr,false);
            }
            else
            {
              mat->getEntry(itn->second, itm->second,cc);
              mat->setEntry(itn->second, itm->second,cc+alpha*complex_t(cmr),false);
            }
          }
        }
        else
        {
          matrix_p->getEntry(i,j,cmc);
          if (cc!=z0)
          {
            mat->getEntry(itn->second, itm->second,cc);
            mat->setEntry(itn->second, itm->second,cc+alpha*cmc,false);
          }
        }
      }
    }
  }
  trace_p->pop();
}

/*! extend storage of matrix when constraints are not local and have column/row combination, assuming scalar matrix entries !
    mat : pointer to the matrix to be reduced
    cdofsr : row dof components
    cdofsc : col dof components
    cu : pointer to the constraints on the unknown u (col)
    cv : pointer to the constraints on test function v (row)
    keepSymmetry : flag to indicate to force symmetry of the extended storage (default not keeping symmetry)
    doRow, doCol, doDiag : flags to indicates which extensions are processed (default all)
*/
void extendStorage(MatrixEntry* mat, std::vector<DofComponent>& cdofsr, std::vector<DofComponent>& cdofsc,
                   const Constraints* cu, const Constraints* cv, bool keepSymmetry, bool doRow, bool doCol, bool doDiag)
{
  trace_p->push("extendStorage(...)");
  if (mat==0 || mat->storageType()==_dense) { trace_p->pop(); return; } //nothing to do !!!
  bool local=true;
  std::map<DofComponent, number_t> melcdofs_u, mrecdofs_u, melcdofs_v, mrecdofs_v;
  //elapsedTime();
  if (cu!=0 && !cu->local)
  {
    cu->extractCdofs(cdofsc, melcdofs_u, mrecdofs_u, false);
    local = (mrecdofs_u.size()==0);
  }
  if (cv!=0 && !cv->local)
  {
    cv->extractCdofs(cdofsr, melcdofs_v, mrecdofs_v, true);
    if (local) local = (mrecdofs_v.size()==0);
  }
  //elapsedTime("extendStorage:extractCdofs");

  if (local) {trace_p->pop(); return;}  //storage extension not required

  //build extended col index
  number_t nbc=mat->nbOfCols(), nbr=mat->nbOfRows();
  std::vector<std::set<number_t> > colindex(nbr);
  std::vector<std::set<number_t> >::iterator itc=colindex.begin();
  MatrixStorage* mst=mat->storagep();
  #ifdef XLIFEPP_WITH_OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (number_t r=0; r<nbr; r++)
    colindex[r] = mst->getCols(r+1); //copy col indices
  //for (number_t r=1; r<=nbr; r++, itc++) *itc = mat->storagep()->getCols(r); //copy col indices
  real_t cr=0.; complex_t cc=0.;
  std::map<DofComponent, number_t>::iterator itm, itn;
  bool storageChange = false;
  //elapsedTime("extendStorage:init colindex");

  const MatrixEntry* cmat=cu->matrixp();
  if (doCol && mrecdofs_u.size()>0)  //column combination
  {
    cmat=cu->matrixp();
    number_t nbe=cu->elcdofs_.size();
    for (itn=melcdofs_u.begin(); itn!=melcdofs_u.end(); itn++)
    {
      std::set<number_t> indrow = mat->storagep()->getRows(itn->second);  //row indices of col
      std::map<DofComponent, number_t>::const_iterator itmd=cu->elcdofs_.find(itn->first);
      number_t i=itmd->second;
      for (itm=mrecdofs_u.begin(); itm!=mrecdofs_u.end(); itm++)
      {
        itmd=cu->recdofs_.find(itm->first);
        number_t j=itmd->second - nbe;
        if (cmat->valueType_ == _real) {cmat->getEntry(i,j,cr); cc=cr;}
        else cmat->getEntry(i,j,cc);
        if (cc!=complex_t(0.)) //add row index in current col
        {
          std::set<number_t>::iterator its;
          for (its=indrow.begin(); its!=indrow.end(); its++) colindex[*its-1].insert(itm->second);
        }
      }
    }
    storageChange=true;
  }
  //elapsedTime("extendStorage: col combination");
  //row combination
  if (doRow)
  {
    cmat=cv->matrixp();
    number_t nbe=cv->elcdofs_.size();
    for (itm=mrecdofs_v.begin(); itm!=mrecdofs_v.end(); itm++)
    {
      std::map<DofComponent, number_t>::const_iterator itmd=cv->recdofs_.find(itm->first);
      if (itmd==cv->recdofs_.end()) itmd=cv->recdofs_.find(itm->first.dual());  //try with dual
      number_t j=itmd->second - nbe;
      for (itn=melcdofs_v.begin(); itn!=melcdofs_v.end(); itn++)
      {
        itmd=cu->elcdofs_.find(itn->first);
        if (itmd==cu->elcdofs_.end()) itmd=cu->elcdofs_.find(itn->first.dual());  //try with dual
        number_t i=itmd->second;
        if (cmat->valueType_ == _real) {cmat->getEntry(i,j,cr); cc=cr;}
        else cmat->getEntry(i,j,cc);
        if (cc!=complex_t(0.)) //add col indices in current col index
        {
          std::set<number_t>& indcol=colindex[itn->second -1];
          colindex[itm->second -1].insert(indcol.begin(),indcol.end());
        }
      }
    }
    storageChange=true;
  }
  //elapsedTime("extendStorage: row combination");
  // insert eliminated x reduced block assuming v_constraints is same as u_constraints
  if (mrecdofs_u.size()>0)
  {
    number_t nbe=melcdofs_u.size();
    for (itn=melcdofs_u.begin(); itn!=melcdofs_u.end(); itn++)
    {
      if (itn->second > colindex.size()) //constraint involve unknowns not present in matrix
        error("ec_unknown_not_found");
      std::map<DofComponent, number_t>::const_iterator itmd=cu->elcdofs_.find(itn->first);
      number_t i=itmd->second;
      for (itm=mrecdofs_u.begin(); itm!=mrecdofs_u.end(); itm++)
      {
        itmd=cu->recdofs_.find(itm->first);
        number_t j=itmd->second - nbe;
        if (cmat->valueType_ == _real)
        {
          cmat->getEntry(i,j,cr);
          if (cr!=0.) colindex[itn->second -1].insert(itm->second);
        }
        else
        {
          cmat->getEntry(i,j,cc);
          if (cc!=0.) colindex[itn->second -1].insert(itm->second);
        }
      }
    }
    storageChange=true;
  }
  //elapsedTime("extendStorage: insert eliminated x reduced block");
  if (storageChange) //create new storage
  {
    std::vector<std::vector<number_t> > vcolindex(nbr);
    std::vector<std::vector<number_t> >::iterator itv=vcolindex.begin();
    itc=colindex.begin();
    for (; itv!=vcolindex.end(); itv++,itc++)
      if (!itc->empty()) *itv = std::vector<number_t>(itc->begin(),itc->end());
    colindex.clear();
    AccessType at = mat->accessType();
    StorageType st= mat->storageType();
    if (at==_sym && !keepSymmetry)
    {
      if (cu != cv || mrecdofs_u.size()>0 || mrecdofs_v.size()>0) at=_dual;
    }
    std::stringstream ss;
    ss << mat->storagep()->stringId << " constraints " << cu << "-" << cv;
    //MatrixStorage* newsto = createMatrixStorage(st,at,nbr,nbc,vcolindex,mat->storagep()->name()+"_extended");
    MatrixStorage* newsto = createMatrixStorage(st,at,nbr,nbc,vcolindex,ss.str());
    //elapsedTime("extendStorage:create storage");
    mat->toStorage(newsto);
    //elapsedTime("extendStorage:move to storage");
    if (at==_dual && !keepSymmetry) mat->symmetry() =_noSymmetry;
  }

  trace_p->pop();
}

/*! full reduction of a pseudo-reduced matrix, assuming scalar matrix entries
    create a new storage from old one by eliminating row en col index related to eliminated dof
    note : cdofsr and cdofsc are not modified !!!

     mat : the matrix to be reduced
     cdofsr : row dof numbering of mat
     cdofsc : col dof numbering of mat
     redcdofsr : row dof numbering of reduced mat
     redcdofsc : col dof numbering of reduced mat
     cu, cv : u and v constraint pointers
*/
void reduceMatrix(MatrixEntry*& mat, std::vector<DofComponent>& cdofsr, std::vector<DofComponent>& cdofsc,
                  std::vector<DofComponent>& redcdofsr, std::vector<DofComponent>& redcdofsc,
                  const Constraints* cu, const Constraints* cv)
{
  trace_p->push("reduceStorage(...)");
  if (mat==0) { trace_p->pop(); return; } //nothing to do !!!
  std::map<DofComponent, number_t> melcdofs_u, mrecdofs_u, melcdofs_v, mrecdofs_v;
  if (cu!=0) cu->extractCdofs(cdofsc, melcdofs_u, mrecdofs_u, false);
  if (cv!=0) cv->extractCdofs(cdofsr, melcdofs_v, mrecdofs_v, true);

  number_t nbe_u = melcdofs_u.size(), nbe_v = melcdofs_v.size();
  number_t nbc=mat->nbOfCols()-nbe_u, nbr=mat->nbOfRows()- nbe_v;
  std::map<DofComponent, number_t>::iterator itm, itn, itme = melcdofs_v.end(), itne=melcdofs_u.end();
  redcdofsr.resize(nbr); redcdofsc.resize(nbc);

  //build map : old col numbering to new col numbering and set redcdofsc/redcdofsr
  number_t i=0;
  std::vector<DofComponent>::iterator itd, itdn=redcdofsc.begin();
  std::vector<number_t> colmap(mat->nbOfCols(),theNumberMax);
  std::vector<number_t> ::iterator ite=colmap.begin();
  for (itd=cdofsc.begin(); itd!=cdofsc.end(); ++itd, ++ite)
    if (melcdofs_u.find(*itd)==itne) {*ite=i++; *itdn++=*itd;}
  itdn=redcdofsr.begin();
  for (itd=cdofsr.begin(); itd!=cdofsr.end(); ++itd)
    if (melcdofs_v.find(*itd)==itme) {*itdn++=*itd;}

  //build reduced col index
  std::vector<std::vector<number_t> > colindex(nbr);
  std::vector<std::vector<number_t> >::iterator itc=colindex.begin();
  std::set<number_t>::iterator its;
  for (number_t r=1; r<=mat->nbOfRows(); r++)
  {
    itm = melcdofs_v.find(cdofsr[r-1]);
    if (itm==itme) // non eliminated row
    {
      std::set<number_t> cols = mat->storagep()->getCols(r);  //get col index
      std::set<number_t> newcols;  //get col index
      for (its=cols.begin(); its!=cols.end(); ++its)
        if ( melcdofs_u.find(cdofsc[*its-1])==itne) newcols.insert(colmap[*its-1]+1);
      itc->insert(itc->begin(),newcols.begin(),newcols.end());
      itc++;
    }
  }
  //create new storage and move matrix to new storage (expansive)
  AccessType at = mat->accessType();
  StorageType sto= mat->storageType();
  if (at==_sym && cu != cv) at=_dual;  //symmetry is lost
  MatrixStorage* newsto = createMatrixStorage(sto,at,nbr,nbc,colindex,mat->storagep()->name()+"_reduced");
  ValueType vt=mat->valueType();
  StrucType st=mat->strucType();
  SymType sy = mat->symmetry();
  MatrixEntry* newmat= new MatrixEntry(vt, st, newsto, mat->nbOfComponents, sy);
  real_t cr; complex_t cc;  i=1;
  for (number_t r=1; r<=mat->nbOfRows(); r++)
  {
    itm = melcdofs_v.find(cdofsr[r-1]);
    if (itm==itme) // non eliminated row
    {
      std::set<number_t> cols = mat->storagep()->getCols(r);  //get col index
      for (its=cols.begin(); its!=cols.end(); ++its)
      {
        number_t j=(*its)-1;
        if ( (r >= *its || at!=_sym || sy==_noSymmetry) && melcdofs_u.find(cdofsc[j])==itne)
        {
          if (vt==_real)   //real matrix
          {
            mat->getEntry(r,*its, cr);
            newmat->setEntry(i,colmap[j]+1,cr);
          }
          else            // complex matrix
          {
            mat->getEntry(r,*its, cc);
            newmat->setEntry(i,colmap[j]+1,cc);
          }
        }
      }
      i++;
    }
  }
  //delete old matrix
  delete mat;
  mat=newmat;
  std::vector<DofComponent> tmp=cdofsc;
  cdofsc=redcdofsc;  redcdofsc=tmp;
  tmp=cdofsr; cdofsr=redcdofsr; redcdofsr=cdofsr;
  trace_p->pop();
}

/*! correct a right hand side b to take into account constraints
    recall that constraints on unknown are of the form  Ue1 + CUr1 = f  (C matrix e1 x r1)
    and constraints on test functions are of the form   Ve2 + DVr2 = 0  (D matrix e2 x r2)
    where ei stands for eliminated indices and ri for reduced indices (e1=e2 and r1=r2 in most cases)
    in the matrix reduction process, a corrector matrix has been computed, say E matrix m x e1
    the correction process consists in
         first step  : b -= E * f_e1       (column combination)
         second step : b_r2 -= Dt * b_e2   (row combination)
         third step  : b_e2 = f_e1         (deletion in case of real reduction)
    in simple cases (Dirichlet for instance), C=D=0 so the second step is not required

    b      : scalar vector to be corrected
    cdofsb : component dof of vector b
    rhsmap : pointer to correction matrix (C)
    cu, cv : pointer to u/v constraints system

*/
void appliedRhsCorrectorTo(VectorEntry* b, const std::vector<DofComponent>& cdofsb,
                           MatrixEntry* rhsmat, const Constraints* cu, const Constraints* cv,
                           const ReductionMethod& rm)
{
  trace_p->push("appliedRhsCorrectorTo()");

  ReductionMethodType meth=rm.method;
  complex_t alpha = rm.alpha;

  if (cu!=0 && cu->rhsp()!=0 && rhsmat!=0) // first step (column combination)
  {
    std::map<DofComponent, number_t>::const_iterator itmd;
    VectorEntry* nrhs;
    const VectorEntry* rhsp=cu->rhsp();
    if (rhsp->valueType_==_real)
    {
      nrhs = new VectorEntry(_real, _scalar,cu->elcdofs_.size());
      number_t k=1; real_t v;
      for (itmd=cu->elcdofs_.begin(); itmd!=cu->elcdofs_.end(); ++itmd, k++)
      {
        rhsp->getEntry(itmd->second,v);
        nrhs->setEntry(k,v);
      }
    }
    else
    {
      nrhs = new VectorEntry(_complex, _scalar,cu->elcdofs_.size());
      number_t k=1; complex_t c;
      for (itmd=cu->elcdofs_.begin(); itmd!=cu->elcdofs_.end(); ++itmd, k++)
      {
        rhsp->getEntry(itmd->second,c);
        nrhs->setEntry(k,c);
      }
    }
    *b -= *rhsmat** nrhs;
    delete nrhs;
  }

  if (cv==0) {trace_p->pop(); return;} //no row reduction

  //create dual cdofsb
  std::vector<DofComponent> dual_cdofsb(cdofsb.size());
  std::vector<DofComponent>::const_iterator itd;
  std::vector<DofComponent>::iterator itdd=dual_cdofsb.begin();
  for (itd=cdofsb.begin(); itd!=cdofsb.end(); itd++, itdd++) *itdd = itd->dual();

  //locate eliminated and reduced dofs in b
  std::map<DofComponent, number_t> belcdofs, brecdofs;
  itdd=dual_cdofsb.begin(); number_t k=1;
  std::map<DofComponent, number_t>::const_iterator itmd;
  for (itd=cdofsb.begin(); itd!=cdofsb.end(); itd++, itdd++, k++)  //goto original row numbering
  {
    itmd=cv->elcdofs_.find(*itd);
    if (itmd!=cv->elcdofs_.end()) belcdofs[itmd->first] = k;
    itmd=cv->elcdofs_.find(*itdd);
    if (itmd!=cv->elcdofs_.end()) belcdofs[itmd->first] = k;
    itmd=cv->recdofs_.find(*itd);
    if (itmd!=cv->recdofs_.end()) brecdofs[itmd->first] = k;
    itmd=cv->recdofs_.find(*itdd);
    if (itmd!=cv->recdofs_.end()) brecdofs[itmd->first] = k;
  }

  if (cv->matrixp()!=0  && cv->matrixp()->nbOfCols()>0)         // second step (row combination)
  {
    number_t nbe=cv->elcdofs_.size();
    VectorEntry rb(b->valueType_,b->strucType_,nbe);

    std::map<DofComponent, number_t>::const_iterator itme;
    for (itmd = belcdofs.begin(); itmd!=belcdofs.end(); itmd++)
    {
      itme=cv->elcdofs_.find(itmd->first);
      if (b->valueType_==_real)
      {
        real_t v;
        b->getEntry(itmd->second,v);
        rb.setEntry(itme->second,v);
      }
      else
      {
        complex_t v;
        b->getEntry(itmd->second, v);
        rb.setEntry(itme->second,v);
      }
    }

    if (brecdofs.size()>0)
    {
      VectorEntry rbe;
      if (cv->matrixp()->valueType_==_complex) rbe = rb * conj(*cv->matrixp());  //warning : constraint matrix is conjugated
      else rbe =  rb** cv->matrixp();

      for (itmd=brecdofs.begin(); itmd!=brecdofs.end(); itmd++)
      {
        itme=cv->recdofs_.find(itmd->first);
        if (b->valueType_==_real)
        {
          real_t v1,v2;
          rbe.getEntry(itme->second-nbe,v1);
          b->getEntry(itmd->second,v2);
          v2-=v1;
          b->setEntry(itmd->second,v2);
        }
        else
        {
          complex_t v1,v2;
          rbe.getEntry(itme->second-nbe,v1);
          b->getEntry(itmd->second,v2);
          v2-=v1;
          b->setEntry(itmd->second,v2);
        }
      }
    }
  }
  if ((meth ==_pseudoReduction || meth==_penalizationReduction)&& cu==cv) // third step : set b eldofs to 0 or u constraints rhs only if cu==cv
    {
      std::map<DofComponent, number_t>::const_iterator itmd2=cv->elcdofs_.begin();
      for (itmd=belcdofs.begin(); itmd!=belcdofs.end(); itmd++, itmd2++)
      {
        if (b->valueType_==_real)
        {
          real_t v=0;
          if (cu->rhsp()!=0) cu->rhsp()->getEntry(itmd2->second,v);
          b->setEntry(itmd->second,alpha.real()*v);
        }
        else
        {
          complex_t v=0;
          if (cu->rhsp()!=0)
          {
            if (cu->rhsp()->valueType_ ==_real)
            {
              real_t r=0;
              cu->rhsp()->getEntry(itmd2->second,r);
              v=r;
            }
            else cu->rhsp()->getEntry(itmd2->second,v);
          }
          b->setEntry(itmd->second,alpha*v);
        }
      }
    }
  else
  {
    if (!(meth ==_pseudoReduction || meth==_penalizationReduction))
      error("reduction_unexpected", words("reduction method",_pseudoReduction), words("reduction method",meth));
  }
  trace_p->pop();
}

// ================================================================================================
/*! apply essential condition to a VectorEntry (user tool)
      v     : the SCALAR vector to be reduced
      cdofs : the component dofs of the vector v
      cs    : the constraints to apply
    once reduced, the constraints reads Ve = f - C*Vr
    where Ve the ne eliminated components, Vr the nr reduced components and C a ne x nr matrix
    this function, first build the vector s = f - C*Vr and then set the eliminated components of v to s
*/
//================================================================================================
void applyEssentialConditions(VectorEntry& v, const std::vector<DofComponent>& cdofs, const Constraints& cs)
{
  trace_p->push("applyEssentialConditions(VectorEntry");
  ValueType vt=v.valueType();
  VectorEntry s;
  if (cs.rhsp()==0)  s=VectorEntry(vt,_scalar,cs.numberOfConstraints());
  else s=*cs.rhsp();
  //deal with constraints matrix if exists
  std::vector<DofComponent>::const_iterator itd;
  std::map<DofComponent, number_t>::const_iterator itrec, itrece;
  const MatrixEntry* matp = cs.matrixp();
  if (matp!=0 && matp->nbOfCols()!=0  && cs.recdofs_.size()>0)
  {
    //locate reduced dofs
    VectorEntry x(vt,_scalar,cs.recdofs_.size());
    number_t i=1, ne=cs.elcdofs_.size();
    itrece=cs.recdofs_.end();
    for (itd=cdofs.begin(); itd!=cdofs.end(); ++itd, i++)
    {
      itrec=cs.recdofs_.find(*itd);
      if (itrec!=itrece) x.setValue(itrec->second-ne,v.getValue(i));
    }
    //do product constraint matrix * x
    s-= *matp * x;
  }
  //locate eliminated dofs and apply constraints to SuTermVector
  number_t i=1;
  itrece=cs.elcdofs_.end();
  for (itd=cdofs.begin(); itd!=cdofs.end(); ++itd, i++)
    {
      itrec=cs.elcdofs_.find(*itd);
      if (itrec!=itrece) v.setValue(i, s.getValue(itrec->second));
    }
  trace_p->pop();
}


// ---------------------------------------------------------------------------------------------------------------------------------
// SetOfConstraints stuff
// ---------------------------------------------------------------------------------------------------------------------------------

// copy, clear, copy constructor, destructor
SetOfConstraints::SetOfConstraints(const SetOfConstraints& soc)
{
  copy(soc);
}

SetOfConstraints& SetOfConstraints::operator = (const SetOfConstraints& soc)
{
  if (this==&soc) return *this;
  clear();
  copy(soc);
  return *this;
}

SetOfConstraints::~SetOfConstraints()
{
  clear();
}

// full copy, has to be cleared before
void SetOfConstraints::copy(const SetOfConstraints& soc)
{
  std::map<const Unknown*, Constraints*>::const_iterator itm=soc.begin();
  for (; itm!=soc.end(); itm++)
  {
    if (itm->second !=0)
    {
      Constraints* c= new Constraints(*itm->second);
      insert(std::make_pair(itm->first,c));
    }
  }
}

//deallocate constraints pointers and reset map to void
void SetOfConstraints::clear()
{
  std::map<const Unknown*, Constraints*>::iterator itm=begin();
  for (; itm!=end(); itm++)
    if (itm->second!=0) delete itm->second;
  std::map<const Unknown*, Constraints*>::clear();
}
// some utilities

Constraints* SetOfConstraints::operator()(const Unknown* u) const
{
  std::map<const Unknown*, Constraints*>::const_iterator itm=find(u);
  if (itm == end()) return 0;
  return itm->second;
}

bool SetOfConstraints::isGlobal() const    //return true if SetOfConstraints has a unique global constraints
{
  if (size()==0) return false;
  if (begin()->first == 0) return true;
  return false;
}

bool SetOfConstraints::isReduced() const   // return true if all constraints in SetOfConstraints have been reduced
{
  std::map<const Unknown*, Constraints*>::const_iterator itm=begin();
  for (; itm!=end(); itm++)
  {
    if (itm->second ==0) return false;
    if (!itm->second->reduced) return false;
  }
  return true;
}

// print utilities
void SetOfConstraints::print(std::ostream& os) const
{
  std::map<const Unknown*, Constraints*>::const_iterator itm=begin();
  for (; itm!=end(); itm++)
  {
    if (itm->first==0) os << "global constraints : ";
    else os << "constraints on unknown " << itm->first->name() << " : ";
    if (itm->second!=0)
    {
      os << eol;
      itm->second->print(os);
    }
    else os << " void !" << eol;
  }
}

std::ostream& operator<<(std::ostream& os, const SetOfConstraints& soc)
{
  soc.print(os);
  return os;
}

} // end of namespace xlifepp

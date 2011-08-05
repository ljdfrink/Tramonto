//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER

#include "dft_TpetraBasicLinProbMgr.hpp"

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_BasicLinProbMgr
(size_t numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm)
  : isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true),
    firstTime_(true),
    numUnknownsPerNode_(numUnknownsPerNode),
    parameterList_(parameterList),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numGlobalBoxNodes_(0),
    numCoarsenedNodes_(0),
    numGlobalCoarsenedNodes_(0),
    comm_(comm),
    curRow_(-1)
{
  return;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_BasicLinProbMgr()
{
  return;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalRowMap
(int *GIDs, const ArrayView<const GlobalOrdinal>& GIDs2, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz)
{
  if (numGlobalNodes_!=0)
  {
    return; // Already been here
  }
  numOwnedNodes_ = GIDs2.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
				  &numOwnedNodes_, &numGlobalNodes_);

  ownedMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs2, 0, comm_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalColMap
(int *GIDs, const ArrayView<const GlobalOrdinal> &GIDs2, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz)
{
  if (numGlobalBoxNodes_!=0)
  {
    return; // Already been here
  }

  numBoxNodes_ = GIDs2.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
			       &numBoxNodes_, &numGlobalBoxNodes_);

  boxMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs2, 0, comm_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setCoarsenedNodesList(const ArrayView<const GlobalOrdinal> &GIDs)
{
  if (numGlobalCoarsenedNodes_!=0)
  {
    return; // Already been here
  }

  numCoarsenedNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
			       &numCoarsenedNodes_, &numGlobalCoarsenedNodes_);

  coarsenedNodesMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs, 0, comm_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeBlockStructure
()
{
  TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Already set block structure.\n");

  const size_t numUnks2 = numOwnedNodes_*numUnknownsPerNode_;
  Array<LocalOrdinal> globalGIDList2(numUnks2);

  // Physics ordering for Basic Linear Problem is natural ordering:
  physicsOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++)
  {
    physicsOrdering_[i] = i;
  }

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++)
  {
    solverOrdering_[physicsOrdering_[i]]=i;
  }

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  ArrayView<const GlobalOrdinal> GIDs2 = ownedMap_->getNodeElementList();
  LocalOrdinal k2=0;
  if (groupByPhysics_)
  {
    for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
    {
      LocalOrdinal ii=physicsOrdering_[i];
      for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
      {
	      globalGIDList2[k2++] = ii*numGlobalNodes_ + GIDs2[j];
      }
    }
  }
  else
  {
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
      {
	LocalOrdinal ii=physicsOrdering_[i];
	globalGIDList2[k2++] = ii + GIDs2[j]*numUnknownsPerNode_;
      }
    }
  }

  globalRowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList2, 0, comm_));
  globalMatrix_ = rcp(new MAT(globalRowMap_, 0));
  globalMatrix_->setObjectLabel("BasicLinProbMgr::globalMatrix");
  globalRhs_ = rcp(new VEC(globalRowMap_));
  globalLhs_ = rcp(new VEC(globalRowMap_));
  ownedToBoxImporter_ = rcp(new IMP(ownedMap_, boxMap_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
  TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error,
		     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error,
		     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");

  isLinearProblemSet_ = false; // We are reinitializing the linear problem

 // AGS: I found that we needed to initialize the matrix even the
 // first time a matrix was filled, because matrix entries are being put
 // in on residual-only fills, which can occur before matrix fills.

  if (!firstTime_)
  {
    globalMatrix_->resumeFill();
    globalMatrix_->setAllToScalar(0.0);
    globalMatrix_->fillComplete();
    globalRhs_->putScalar(0.0);
    globalLhs_->putScalar(0.0);
  }
 }
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRhsValue
(GlobalOrdinal ownedPhysicsID, GlobalOrdinal ownedNode, Scalar value)
{
  LocalOrdinal rhsLID = ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  globalRhs_->sumIntoLocalValue(rhsLID, value);
  //cout << std::setprecision(2);
  //cout << "b[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode<<"] = " << value << endl;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, Scalar value)
{
  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);

  //cout << std::setprecision(2);
  //cout << "A[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode
  //     << "][boxPhysicsID="  <<boxPhysicsID  <<"][boxNode="  <<boxNode
  //     << "][rowGID="        <<rowGID        <<"][colGID="   <<colGID
  //     << "] = " << value << endl;
  if (firstTime_) {
    if (rowGID!=curRow_) {
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else {
    Array<Scalar> vals(1);
    vals[0] = value;
    Array<GlobalOrdinal> globalCols(1);
    globalCols[0] = colGID;
    globalMatrix_->sumIntoGlobalValues(rowGID, globalCols, vals);
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  if (curRowValues_.empty()) return;

  typename std::map<LocalOrdinal, Scalar>::iterator pos2;
  for (pos2 = curRowValues_.begin(); pos2 != curRowValues_.end(); ++pos2) {
    indices_.append(pos2->first);
    values_.append(pos2->second);
  }

  globalMatrix_->insertGlobalValues(curRow_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode)
{
  TEST_FOR_EXCEPTION(globalMatrix_.get()==0, std::runtime_error, "Global Matrix is not constructed, must set debug flag to enable this feature.\n");

  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);
  size_t numEntries;
  ArrayView<const GlobalOrdinal> indices;
  ArrayView<const Scalar> values;
  if (globalMatrix_->isGloballyIndexed())
  {
    globalMatrix_->getGlobalRowView(rowGID, indices, values); // get view of current row
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; i++)
    {
      if (colGID==indices[i])
      {
	return(values[i]);
      }
    }
  }
  else
  {
    rowGID = globalMatrix_->getRowMap()->getLocalElement(rowGID); // get local row ID
    colGID = globalMatrix_->getColMap()->getLocalElement(colGID); // get local column ID
    globalMatrix_->getLocalRowView(rowGID, indices, values); // get view of current row
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; i++)
    {
      if (colGID==indices[i])
      {
	return(values[i]);
      }
    }
  }

  return(0.0);
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID,
 const ArrayView<const LocalOrdinal>& boxNodeList, const ArrayView<const Scalar>& values)
{
  for (LocalOrdinal i=0; i<boxNodeList.size(); i++)
  {
    insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNodeList[i], values[i]);
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, const ArrayView<const LocalOrdinal> &boxPhysicsIDList,
 LocalOrdinal boxNode, const ArrayView<const Scalar> &values)
{
  for (LocalOrdinal i=0; i<boxPhysicsIDList.size(); i++)
  {
    insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsIDList[i], boxNode, values[i]);
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
    return; // nothing to do

  if (firstTime_) {
    insertRow();

    if(!globalMatrix_->isFillComplete()){
      globalMatrix_->fillComplete();
    }
  }

  isLinearProblemSet_ = true;
  firstTime_ = false;

  //writeMatrix("basica.dat", "", "");
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setRhs
(const ArrayView<const ArrayView<const Scalar> >& b)
{
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      globalRhs_->replaceLocalValue(ownedToSolverLID(i,j), b[i][j]);
    }
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setLhs
(const ArrayView<const ArrayView<const Scalar> > &x) const
{
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    ArrayRCP<const Scalar> xtmp = exportC2R(Teuchos::arcpFromArrayView<const Scalar>(x[i])); // Use simple import
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      globalLhs_->replaceLocalValue(ownedToSolverLID(i,j), xtmp[j]);
    }
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getLhs
() const
{
  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalLhs_->get1dViewNonConst();

  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    ArrayRCP<Scalar> temp(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      temp[j] = tmp[ownedToSolverLID(i,j)];
    }
    ArrayOfPtrs[i] = importR2C(temp.getConst()); // Use simple import (uses view so doesn't work)
  }
  return ArrayOfPtrs;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getRhs
() const
{
  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalRhs_->get1dViewNonConst();
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    ArrayOfPtrs[i] = Teuchos::arcp<Scalar>(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      ArrayOfPtrs[i][j] = tmp[ownedToSolverLID(i,j)];
    }
  }
  return ArrayOfPtrs;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setupSolver
()
{
  TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
		     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");

#ifdef SUPPORTS_STRATIMIKOS
  thyraRhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalRhs_);
  thyraLhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalLhs_);
  thyraOp_ = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalMatrix_);

  solver_ = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver_->readParameters(out.get());

  lowsFactory_ = solver_->createLinearSolveStrategy("");
  lows_ = linearOpWithSolve<Scalar>(*lowsFactory_, thyraOp_);
#else
  //printf("Node num rows = %d\n", globalMatrix_->getNodeNumRows());
  //printf("Node num cols = %d\n", globalMatrix_->getNodeNumCols());
  problem_ = rcp(new LinPROB(globalMatrix_, globalLhs_, globalRhs_));
  // No preconditioner for now
  // Ifpack2::Factory factory;
  // RCP<const MAT> const_globalMatrix_ = Teuchos::rcp_implicit_cast<const MAT>(globalMatrix_);
  // preconditioner_ = factory.create("ILUT", const_globalMatrix_);
  // preconditioner_->setParameters(*parameterList_);
  // preconditioner_->initialize();
  // preconditioner_->compute();
  // problem_->setLeftPrec(preconditioner_);
  TEST_FOR_EXCEPT(problem_->setProblem() == false);
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, parameterList_));
#endif
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
solve
()
{
#ifdef KDEBUG
  printf("\n\n\n\ndft_BasicLinProbMgr::solve()\n\n\n\n");
#endif

#ifdef SUPPORTS_STRATIMIKOS
  SolveStatus<Scalar> status = lows_->solve(Thyra::NOTRANS, *thyraRhs_, thyraLhs_.ptr());
#else
  ReturnType ret = solver_->solve();
#endif

  bool writeMatrixNow = false;
  if (writeMatrixNow)
  {
    writeMatrix("A_ref.dat", "GlobalMatrix", "GlobalMatrix");
    writeLhs("x_ref.dat");
    writeRhs("b_ref.dat");
    writePermutation("p_ref.dat");
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const
{
  setLhs(x);
  globalMatrix_->apply(*globalLhs_, *globalRhs_);

  return(getRhs());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const ArrayRCP<const ArrayRCP<const Scalar> >& xOwned) const
{
  ArrayRCP<ArrayRCP<Scalar> > my_xBox = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    my_xBox[i] = importR2C(xOwned[i]);
  }

  return(my_xBox);
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const ArrayRCP<const Scalar> &aOwned) const
{
  RCP<VEC> owned2 = rcp(new VEC(ownedMap_, aOwned()));
  RCP<VEC>  box2 = rcp(new VEC(boxMap_));

  box2->doImport(*owned2, *ownedToBoxImporter_, Tpetra::INSERT);

  return (box2->get1dViewNonConst());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
exportC2R
(const ArrayRCP<const Scalar>& aBox) const
{
  RCP<VEC> owned2 =  rcp(new VEC(ownedMap_));
  RCP<VEC> box = rcp(new VEC(boxMap_, aBox()));

  owned2->doExport(*box, *ownedToBoxImporter_, Tpetra::INSERT); // Use importer, but zero out off-processor contributions.

  return(owned2->get1dViewNonConst());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeMatrix
(const char * filename, const char * matrixName, const char * matrixDescription) const  {
  //int dft_BasicLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
  return;
      //(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_, matrixName, matrixDescription));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeLhs
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writeLhs(const char * filename) const  {
  return;
      //(EpetraExt::MultiVectorToMatlabFile(filename, *globalLhs_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeRhs
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writeRhs(const char * filename) const  {
  return;
      //(EpetraExt::MultiVectorToMatlabFile(filename, *globalRhs_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writePermutation
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writePermutation(const char * filename) const  {
  return;
    //(EpetraExt::BlockMapToMatrixMarketFile(filename, *globalRowMap_, " ", " ", false));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
checkPhysicsOrdering()
const  {
  TEST_FOR_EXCEPTION(physicsOrdering_.size()==0, std::runtime_error, "No unknowns are registered with this problem manager.\n");

  size_t numUnks = physicsOrdering_.size();
  Array<Scalar> tmp(numUnks);
  for (LocalOrdinal i=0; i<numUnks; i++)
  {
    LocalOrdinal curID = physicsOrdering_[i];
    TEST_FOR_EXCEPTION(curID <0, std::runtime_error, "Invalid unknown number " << curID << " is less than 0.\n");
    TEST_FOR_EXCEPTION(curID>=numUnks, std::runtime_error, "Invalid unknown number " << curID << " is greater than or equal to the number of unknowns (" << numUnks << ").\n");
      tmp[curID] = tmp[curID]+1;
      // Increment counter for this ID (at the end each ID should appear exactly one time).
  }

  for (LocalOrdinal i=0; i<numUnks; i++)
  {
    TEST_FOR_EXCEPTION(tmp[i]==0, std::runtime_error, "Unknown number " << i << " is not present and should be.\n");
    TEST_FOR_EXCEPTION(tmp[i]>1, std::runtime_error, "Unknown number " << i << " is present " << tmp[i] << " times and should be present only once.\n");
  }

}
template class dft_BasicLinProbMgr<double, int, int>;

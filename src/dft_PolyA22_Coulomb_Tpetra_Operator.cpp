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

#include "dft_PolyA22_Coulomb_Tpetra_Operator.hpp"

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA22_Coulomb_Tpetra_Operator
(const RCP<const MAP> & cmsMap, const RCP<const MAP> & densityMap,
 const RCP<const MAP> & poissonMap, const RCP<const MAP> & cmsDensMap,
 const RCP<const MAP> & block2Map,  RCP<ParameterList> parameterList)
  : dft_PolyA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>
	      (cmsMap, densityMap, cmsDensMap, parameterList),
    poissonMap_(poissonMap),
    cmsDensMap_(cmsDensMap),
    block2Map_(block2Map),
    curPoissonRow_(-1),
    curCPRow_(-1),
    curPDRow_(-1)
{
  poissonOnPoissonMatrix_ = rcp(new MAT(poissonMap, 0));
  poissonOnPoissonMatrixOp_ = rcp(new MMOP(poissonOnPoissonMatrix_));
  cmsOnPoissonMatrix_ = rcp(new MAT(cmsMap, 0));
  cmsOnPoissonMatrixOp_ = rcp(new MMOP(cmsOnPoissonMatrix_));
  poissonOnDensityMatrix_ = rcp(new MAT(poissonMap, 0));
  poissonOnDensityMatrixOp_ = rcp(new MMOP(poissonOnDensityMatrix_));
  Label_ = "dft_PolyA22_Coulomb_Tpetra_Operator";
  cmsOnDensityMatrix_->setObjectLabel("PolyA22Coulomb::cmsOnDensityMatrix");
  cmsOnCmsMatrix_->setObjectLabel("PolyA22Coulomb::cmsOnCmsMatrix");
  poissonOnPoissonMatrix_->setObjectLabel("PolyA22Coulomb::poissonOnPoissonMatrix");
  cmsOnPoissonMatrix_->setObjectLabel("PolyA22Coulomb::cmsOnPoissonMatrix");
  poissonOnDensityMatrix_->setObjectLabel("PolyA22Coulomb::poissonOnDensityMatrix");

} //end constructor
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA22_Coulomb_Tpetra_Operator
()
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{

  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
  {
    if (!isFLinear_)
    {
      cmsOnDensityMatrix_->resumeFill();
      cmsOnDensityMatrix_->setAllToScalar(STMS::zero());
    } //end if
    cmsOnCmsMatrix_->resumeFill();
    cmsOnCmsMatrix_->setAllToScalar(STMS::zero());
    densityOnDensityMatrix_->putScalar(STS::zero());
    densityOnDensityInverse_->putScalar(STS::zero());
    densityOnCmsMatrix_->resumeFill();
    densityOnCmsMatrix_->setAllToScalar(STMS::zero());
    poissonOnPoissonMatrix_->resumeFill();
    poissonOnPoissonMatrix_->setAllToScalar(STMS::zero());
    cmsOnPoissonMatrix_->resumeFill();
    cmsOnPoissonMatrix_->setAllToScalar(STMS::zero());
    poissonOnDensityMatrix_->resumeFill();
    poissonOnDensityMatrix_->setAllToScalar(STMS::zero());
  } //end if
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, MatScalar value, GlobalOrdinal blockColFlag)
{
  // if poisson then blockColFlag = 0
  // if density then blockColFlag = 1
  // if cms then blockColFlag = 2

  Array<GlobalOrdinal> cols(1);
  cols[0] = colGID;
  Array<MatScalar> vals(1);
  vals[0] = value;

  /* The poissonMatrix_, poissonOnDensityMatrix_, cmsOnPoissonMatrix_, and cmsOnDensityMatrix_ values do not change between iterations */

  if (poissonMap_->isNodeGlobalElement(rowGID)) { // Insert into poissonOnPoissonMatrix or poissonOnDensityMatrix
    if ( blockColFlag == 0 ) { // Insert into poissonOnPoissonMatrix
      if (firstTime_) {
	if (rowGID != curRow_) {
	  insertRow();
	  curRow_ = rowGID;
	}
	curRowValuesPoissonOnPoisson_[colGID] += value;
      }
      else
	poissonOnPoissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    }
    else if ( blockColFlag == 1 ) { // Insert into poissonOnDensityMatrix
      if (firstTime_) {
      	if (rowGID != curRow_) {
      	  insertRow();
      	  curRow_ = rowGID;
      	}
      	curRowValuesPoissonOnDensity_[colGID] += value;
      }
      else
      	poissonOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- row in poissonMap, but blockColFlag not set for Poisson or density equations.");
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } //end if poissonMap_.MyGID(rowGID)
  else if (cmsMap_->isNodeGlobalElement(rowGID)) { // Insert into cmsOnPoissonMatrix or cmsOnCmsMatrix or cmsOnDensityMatrix
    if ( blockColFlag == 0 ) { // Insert into cmsOnPoissonMatrix
      if (firstTime_) {
	if (rowGID != curRow_) {
	  insertRow();
	  curRow_ = rowGID;
	}
	curRowValuesCmsOnPoisson_[colGID] += value;
      }
      else
	cmsOnPoissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    }
    else if ( blockColFlag == 2 ) { // Insert into cmsOnCmsMatrix
      if (firstTime_) {
	if (rowGID!=curRow_) {
	  insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValuesCmsOnCms_[colGID] += value;
      }
      else
	cmsOnCmsMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    }
    else if (blockColFlag == 1) { // Insert into cmsOnDensityMatrix ("F matrix")
      if (firstTime_) {
      	if (rowGID!=curRow_) {
      	  insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      	  curRow_=rowGID;
      	}
      	curRowValuesCmsOnDensity_[colGID] += value;
      }
      else if (!isFLinear_) {
      	cmsOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      }
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- row in cmsMap, but blockColFlag not set for Poisson,density, or cms equations.");
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into cmsOnPoisson or cmsOnCmsMatrix or cmsOnDensityMatrix
  else if (densityMap_->isNodeGlobalElement(rowGID)) { // Insert into densityOnDensityMatrix or densityOnCmsMatrix
    if ( blockColFlag == 1 ) { // Insert into densityOnDensityMatrix
      if (rowGID!=colGID) {
	char err_msg[200];
	sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- Inserting non-diagonal element into densityOnDensity matrix.");
	TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg); // Confirm that this is a diagonal value
      }
      densityOnDensityMatrix_->sumIntoLocalValue(densityMap_->getLocalElement(rowGID), value);
    }
    else if ( blockColFlag == 2) { // Insert into densityOnCmsMatrix
      // The density-on-cms matrix is diagonal for most but not all use cases.
      if (firstTime_) {
	if (rowGID!=curRow_) {
	  insertRow();  // Dump the current contents of curRowValues maps  into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValuesDensityOnCms_[colGID] += value;
      }
      else
      	densityOnCmsMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- row in densityMap, but blockColFlag not set for cms or density equations.");
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into densityOnDensityMatrix or densityOnCmsMatrix
  else { // Problem! rowGID not in cmsMap or densityMap or poissonMap
    char err_msg[200];
    sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): rowGID=%i not in cmsMap,densityMap, or poissonMap.",rowGID);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
  }

} //end insertMatrixValue
//=============================================================================
template<class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  // Fill matrix rows
  if (!curRowValuesCmsOnDensity_.empty()) {
    int numEntriesCmsOnDensity = curRowValuesCmsOnDensity_.size();
    if (numEntriesCmsOnDensity>indicesCmsOnDensity_.size()) {
      indicesCmsOnDensity_.resize(numEntriesCmsOnDensity);
      valuesCmsOnDensity_.resize(numEntriesCmsOnDensity);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesCmsOnDensity_.begin(); pos != curRowValuesCmsOnDensity_.end(); ++pos) {
      indicesCmsOnDensity_[i] = pos->first;
      valuesCmsOnDensity_[i++] = pos->second;
    }
    cmsOnDensityMatrix_->insertGlobalValues(curRow_, indicesCmsOnDensity_, valuesCmsOnDensity_);
  }
  if (!curRowValuesCmsOnCms_.empty()) {
    int numEntriesCmsOnCms = curRowValuesCmsOnCms_.size();
    if (numEntriesCmsOnCms>indicesCmsOnCms_.size()) {
      indicesCmsOnCms_.resize(numEntriesCmsOnCms);
      valuesCmsOnCms_.resize(numEntriesCmsOnCms);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesCmsOnCms_.begin(); pos != curRowValuesCmsOnCms_.end(); ++pos) {
      indicesCmsOnCms_[i] = pos->first;
      valuesCmsOnCms_[i++] = pos->second;
    }
    cmsOnCmsMatrix_->insertGlobalValues(curRow_, indicesCmsOnCms_, valuesCmsOnCms_);
  }
  if (!curRowValuesPoissonOnPoisson_.empty()) {
    int numEntriesPoissonOnPoisson = curRowValuesPoissonOnPoisson_.size();
    if (numEntriesPoissonOnPoisson>indicesPoissonOnPoisson_.size()) {
      indicesPoissonOnPoisson_.resize(numEntriesPoissonOnPoisson);
      valuesPoissonOnPoisson_.resize(numEntriesPoissonOnPoisson);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesPoissonOnPoisson_.begin(); pos != curRowValuesPoissonOnPoisson_.end(); ++pos) {
      indicesPoissonOnPoisson_[i] = pos->first;
      valuesPoissonOnPoisson_[i++] = pos->second;
    }
    poissonOnPoissonMatrix_->insertGlobalValues(curRow_, indicesPoissonOnPoisson_, valuesPoissonOnPoisson_);
  }
  if (!curRowValuesCmsOnPoisson_.empty()) {
    int numEntriesCmsOnPoisson = curRowValuesCmsOnPoisson_.size();
    if (numEntriesCmsOnPoisson>indicesCmsOnPoisson_.size()) {
      indicesCmsOnPoisson_.resize(numEntriesCmsOnPoisson);
      valuesCmsOnPoisson_.resize(numEntriesCmsOnPoisson);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesCmsOnPoisson_.begin(); pos != curRowValuesCmsOnPoisson_.end(); ++pos) {
      indicesCmsOnPoisson_[i] = pos->first;
      valuesCmsOnPoisson_[i++] = pos->second;
    }
    cmsOnPoissonMatrix_->insertGlobalValues(curRow_, indicesCmsOnPoisson_, valuesCmsOnPoisson_);
  }
  if (!curRowValuesPoissonOnDensity_.empty()) {
    int numEntriesPoissonOnDensity = curRowValuesPoissonOnDensity_.size();
    if (numEntriesPoissonOnDensity>indicesPoissonOnDensity_.size()) {
      indicesPoissonOnDensity_.resize(numEntriesPoissonOnDensity);
      valuesPoissonOnDensity_.resize(numEntriesPoissonOnDensity);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesPoissonOnDensity_.begin(); pos != curRowValuesPoissonOnDensity_.end(); ++pos) {
      indicesPoissonOnDensity_[i] = pos->first;
      valuesPoissonOnDensity_[i++] = pos->second;
    }
    poissonOnDensityMatrix_->insertGlobalValues(curRow_, indicesPoissonOnDensity_, valuesPoissonOnDensity_);
  }

  if (!curRowValuesDensityOnCms_.empty()) {
    size_t numEntriesDensityOnCms = curRowValuesDensityOnCms_.size();
    if (numEntriesDensityOnCms>indicesDensityOnCms_.size()) {
      indicesDensityOnCms_.resize(numEntriesDensityOnCms);
      valuesDensityOnCms_.resize(numEntriesDensityOnCms);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesDensityOnCms_.begin(); pos != curRowValuesDensityOnCms_.end(); ++pos) {
      indicesDensityOnCms_[i] = pos->first;
      valuesDensityOnCms_[i++] = pos->second;
    }
    densityOnCmsMatrix_->insertGlobalValues(curRow_, indicesDensityOnCms_, valuesDensityOnCms_);
  }

  indicesCmsOnDensity_.clear();
  valuesCmsOnDensity_.clear();
  indicesCmsOnCms_.clear();
  valuesCmsOnCms_.clear();
  indicesPoissonOnPoisson_.clear();
  valuesPoissonOnPoisson_.clear();
  indicesCmsOnPoisson_.clear();
  valuesCmsOnPoisson_.clear();
  indicesPoissonOnDensity_.clear();
  valuesPoissonOnDensity_.clear();
  indicesDensityOnCms_.clear();
  valuesDensityOnCms_.clear();
  curRowValuesCmsOnDensity_.clear();
  curRowValuesCmsOnCms_.clear();
  curRowValuesPoissonOnPoisson_.clear();
  curRowValuesCmsOnPoisson_.clear();
  curRowValuesPoissonOnDensity_.clear();
  curRowValuesDensityOnCms_.clear();

}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
  {
    return;
  }
  insertRow(); // Dump any remaining entries

  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_->getGlobalNumElements()==densityMap_->getGlobalNumElements();

  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
  cmsOnCmsMatrix_->fillComplete(pl);

  if (!isFLinear_) {
    cmsOnDensityMatrix_->fillComplete(densityMap_, cmsMap_, pl);
  }

  poissonOnPoissonMatrix_->fillComplete(pl);
  cmsOnPoissonMatrix_->fillComplete(poissonMap_, cmsMap_, pl);
  poissonOnDensityMatrix_->fillComplete(densityMap_, poissonMap_, pl);

  if (!hasDensityOnCms)  // Confirm that densityOnCmsMatrix is zero
  {
    //    Scalar normvalue = densityOnCmsMatrix_->normInf();
    //    TEUCHOS_TEST_FOR_EXCEPT(normvalue!=STS::zero());
  } else {
    insertRow(); // Dump any remaining entries
    densityOnCmsMatrix_->fillComplete(cmsMap_, densityMap_, pl);
  }

  // Form the inverse of the densityOnDensityMatrix
  densityOnDensityInverse_->reciprocal(*densityOnDensityMatrix_);

  // Use a diagonal preconditioner for the cmsOnCmsMatrix
  RCP<const MAT> const_matrix = Teuchos::rcp_implicit_cast<const MAT>(cmsOnCmsMatrix_);
  cmsOnCmsInverse_ = rcp(new DIAGONAL(const_matrix));
  cmsOnCmsInverseOp_ = rcp(new DIAGONAL_OP(cmsOnCmsInverse_));
  TEUCHOS_TEST_FOR_EXCEPT(cmsOnCmsInverse_==Teuchos::null);
  cmsOnCmsInverse_->initialize();
  cmsOnCmsInverse_->compute();

  if (firstTime_) {
#if ENABLE_MUELU == 1
#if LINSOLVE_PREC_DOUBLE_DOUBLE == 1 || LINSOLVE_PREC_QUAD_DOUBLE == 1
    // Default of SuperLU doesn't compile with double-double or quad-double, so use ILUT instead.
    mueluPP_ = rcp(new Xpetra::TpetraCrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node, typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >(poissonOnPoissonMatrix_));
    mueluPP  = rcp(new Xpetra::CrsMatrixWrap<MatScalar, LocalOrdinal, GlobalOrdinal, Node, typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>(mueluPP_));
    mueluList_ = parameterList_->get("mueluList");

    MueLu::MLParameterListInterpreter<MatScalar, LocalOrdinal, GlobalOrdinal, Node> mueluFactory(mueluList_);

    H_ = mueluFactory.CreateHierarchy();
    H_->setVerbLevel(Teuchos::VERB_NONE);
    H_->GetLevel(0)->Set("A", mueluPP);
    ParameterList coarseParamList;
    coarseParamList.set("fact: level-of-fill", 0);
    RCP<SmootherPrototype> coarsePrototype = rcp( new TrilinosSmoother("ILUT", coarseParamList) );
    RCP<SmootherFactory> coarseSolverFact = rcp( new SmootherFactory(coarsePrototype, Teuchos::null) );
    RCP<FactoryManager> fm = rcp( new FactoryManager() );
    fm->SetFactory("CoarseSolver", coarseSolverFact);
    H_->Setup(*fm);
#else
    // Okay to use default of SuperLU if not double-double or quad-double
    mueluPP_ = rcp(new Xpetra::TpetraCrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node, typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >(poissonOnPoissonMatrix_));
    mueluPP = rcp(new Xpetra::CrsMatrixWrap<MatScalar, LocalOrdinal, GlobalOrdinal, Node, typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>(mueluPP_));
    mueluList_ = parameterList_->sublist("mueluList");

    MueLu::MLParameterListInterpreter<MatScalar, LocalOrdinal, GlobalOrdinal, Node> mueluFactory(mueluList_);

    H_ = mueluFactory.CreateHierarchy();
    H_->setVerbLevel(Teuchos::VERB_NONE);
    H_->GetLevel(0)->Set("A", mueluPP);
    mueluFactory.SetupHierarchy(*H_);
#endif
    poissonOnPoissonInverse_ = rcp(new MueLu::TpetraOperator<MatScalar, LocalOrdinal, GlobalOrdinal, Node>(H_));
    poissonOnPoissonInverseMixed_ = rcp(new MOP((RCP<OP_M>)poissonOnPoissonInverse_));
#endif
  }

  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{
  // If F is in SW (F_location_ == 0):
  // The true A22 block is of the form:

  // |  P       Pd      0    |
  // |  0       Ddd     Ddc  |
  // |  Cp      F       Dcc  |

  // where
  // P is Poisson matrix,
  // Pd is Poisson on Density matrix,
  // Ddd is Density on Density (diagonal),
  // Ddc is Density on Cms (diagonal with small coefficient values),
  // Cp is Cms on Poisson matrix,
  // F is Cms on Density (fairly dense),
  // Dcc is Cms on Cms (diagonal).
  //
  // We will approximate A22 with:

  // |  P       Pd      0    |
  // |  0       Ddd     0    |
  // |  Cp      F       Dcc  |

  // replacing Ddc with a zero matrix for the applyInverse method only.

  // Our algorithm is then:
  // Y1 = Ddd \ X1
  // Y0 = P \ (X0 - Pd*Y1)
  // Y2 = Dcc \ (X2 - Cp*Y0 - F*Y1)

  // A similar algorithm is found when F is in the NE quadrant

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  printf("\n\n\n\ndft_PolyA22_Coulomb_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements(); // == numCmsElements
  size_t numPoissonElements = poissonMap_->getNodeNumElements();

  // X0 is a view of the first numPoisson elements of X
  RCP<const MV> X0 = X.offsetView(poissonMap_, 0);
  // X1 is a view of the middle numDensity/numCms elements of X
  RCP<const MV> X1;
  // X2 is a view of the last numDensity/numCms elements of X
  RCP<const MV> X2;
  // Y0 is a view of the first numPoisson elements of Y
  RCP<MV> Y0 = Y.offsetViewNonConst(poissonMap_, 0);
  // Y1 is a view of the middle numDensity/numCms elements of Y
  RCP<MV> Y1;
  // Y2 is a view of the last numDensity/numCms elements of Y
  RCP<MV> Y2;

  Scalar ONE = STS::one();
  Scalar ZERO = STS::zero();

  if (F_location_ == 1)
  {
    //F in NE part
    X1 = X.offsetView(cmsMap_, numPoissonElements);
    Y1 = Y.offsetViewNonConst(cmsMap_, numPoissonElements);
    X2 = X.offsetView(densityMap_, numPoissonElements+numCmsElements);
    Y2 = Y.offsetViewNonConst(densityMap_, numPoissonElements+numCmsElements);
  }
  else
  {
    //F in SW part
    X1 = X.offsetView(densityMap_, numPoissonElements);
    Y1 = Y.offsetViewNonConst(densityMap_, numPoissonElements);
    X2 = X.offsetView(cmsMap_, numPoissonElements+numDensityElements);
    Y2 = Y.offsetViewNonConst(cmsMap_, numPoissonElements+numDensityElements);
  }

  RCP<MV > Y0tmp = rcp(new MV(*Y0));
  RCP<MV > Y1tmp1 = rcp(new MV(*Y1));
  RCP<MV > Y1tmp2 = rcp(new MV(*Y1));
  RCP<MV > Y2tmp1 = rcp(new MV(*Y2));
  RCP<MV > Y2tmp2 = rcp(new MV(*Y2));

  if (F_location_ == 1)
  {
    // Third block row: Y2 = DD\X2
    Y2->elementWiseMultiply(ONE, *densityOnDensityInverse_, *X2, ZERO);

    // First block row: Y0 = PP \ (X0 - PD*Y2);
    poissonOnDensityMatrixOp_->apply(*Y2, *Y0tmp);
    Y0tmp->update(ONE, *X0, -ONE);
    Y0->putScalar(ZERO);
#if ENABLE_MUELU == 1
    poissonOnPoissonInverseMixed_->apply(*Y0tmp, *Y0);
#endif

    // Second block row: Y1 = CC \ (X1 - CP*Y0 - CD*Y2)
    cmsOnPoissonMatrixOp_->apply(*Y0, *Y1tmp1);
    cmsOnDensityMatrixOp_->apply(*Y2, *Y1tmp2);
    Y1tmp1->update(ONE, *X1, -ONE, *Y1tmp2, -ONE);
    cmsOnCmsInverseOp_->apply(*Y1tmp1, *Y1);

  }
  else
  {
    // Second block row: Y1 = DD\X1
    Y1->elementWiseMultiply(ONE, *densityOnDensityInverse_, *X1, ZERO);

    // First block row: Y0 = PP \ (X0 - PD*Y1);
    poissonOnDensityMatrixOp_->apply(*Y1, *Y0tmp);
    Y0tmp->update(ONE, *X0, -ONE);
    Y0->putScalar(ZERO);
#if ENABLE_MUELU == 1
    poissonOnPoissonInverseMixed_->apply(*Y0tmp, *Y0);
#endif

    // Third block row: Y2 = CC \ (X2 - CP*Y0 - CD*Y1)
    cmsOnPoissonMatrixOp_->apply(*Y0, *Y2tmp1);
    cmsOnDensityMatrixOp_->apply(*Y1, *Y2tmp2);
    Y2tmp1->update(ONE, *X2, -ONE, *Y2tmp2, -ONE);
    cmsOnCmsInverseOp_->apply(*Y2tmp1, *Y2);

  }
  /*
#ifdef SUPPORTS_STRATIMIKOS
  RCP<ThyraMV> thyraY = createMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Y0);
  RCP<ThyraOP> thyraOp = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(poissonMatrix_);

  RCP<DefaultLinearSolverBuilder> solver = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver->readParameters(out.get());

  RCP<ThyraLOWSFactory> lowsFactory = solver->createLinearSolveStrategy("");
  RCP<ThyraLOWS> lows = linearOpWithSolve<Scalar>(*lowsFactory, thyraOp);
  SolveStatus<Scalar> status = lows->solve(Thyra::NOTRANS, *thyraY, thyraY.ptr());
#else
  RCP<LinPROB> problem = rcp(new LinPROB(poissonMatrix_, Y0, Y0));
  TEUCHOS_TEST_FOR_EXCEPT(problem->setProblem() == false);
  RCP<SolMGR> solver = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem, parameterList_));
  ReturnType ret = solver->solve();
#endif
  */

} //end applyInverse
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements();
  size_t numPoissonElements = poissonMap_->getNodeNumElements();

  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_->getGlobalNumElements()==densityMap_->getGlobalNumElements();

  // X0 is a view of the first numPoisson elements of X
  RCP<const MV> X0 = X.offsetView(poissonMap_, 0);
  // X1 is a view of the middle numDensityElements/numCms elements of X
  RCP<const MV> X1;
  // X2 is a view of the last numDensity/numCms elements of X
  RCP<const MV> X2;
  // Y0 is a view of the first numPoisson elements of Y
  RCP<MV > Y0 = Y.offsetViewNonConst(poissonMap_, 0);
  // Y1 is a view of the middle numDensity/numCms elements of Y
  RCP<MV > Y1;
  // Y2 is a view of the last numDensity/numCms elements of Y
  RCP<MV > Y2;

  Scalar ONE = STS::one();
  Scalar ZERO = STS::zero();

  if (F_location_ == 1)
  {
    // F in NE
    X1 = X.offsetView(cmsMap_, numPoissonElements);
    Y1 = Y.offsetViewNonConst(cmsMap_, numPoissonElements);
    X2 = X.offsetView(densityMap_, numPoissonElements+numCmsElements);
    Y2 = Y.offsetViewNonConst(densityMap_, numPoissonElements+numCmsElements);
  }
  else
  {
    // F in SW
    X1 = X.offsetView(densityMap_, numPoissonElements);
    Y1 = Y.offsetViewNonConst(densityMap_, numPoissonElements);
    X2 = X.offsetView(cmsMap_, numPoissonElements+numDensityElements);
    Y2 = Y.offsetViewNonConst(cmsMap_, numPoissonElements+numDensityElements);
  }

  RCP<MV > Y0tmp = rcp(new MV(*Y0));
  RCP<MV > Y0tmp2 = rcp(new MV(*Y0));
  RCP<MV > Y1tmp1 = rcp(new MV(*Y1));
  RCP<MV > Y1tmp2 = rcp(new MV(*Y1));
  RCP<MV > Y2tmp1 = rcp(new MV(*Y2));
  RCP<MV > Y2tmp2 = rcp(new MV(*Y2));

  if (F_location_ == 1)
  {
    // First block row
    poissonOnPoissonMatrixOp_->apply(*X0, *Y0);
    poissonOnDensityMatrixOp_->apply(*X2, *Y0tmp);
    Y0->update(ONE, *Y0tmp, ONE);

    // Second block row
    cmsOnPoissonMatrixOp_->apply(*X0, *Y1);
    cmsOnCmsMatrixOp_->apply(*X1, *Y1tmp1);
    cmsOnDensityMatrixOp_->apply(*X2, *Y1tmp2);
    Y1->update(ONE, *Y1tmp1, ONE);
    Y1->update(ONE, *Y1tmp2, ONE);

    // Third block row
    if (hasDensityOnCms) {
      densityOnCmsMatrixOp_->apply(*X1, *Y2);
      Y2->elementWiseMultiply(ONE, *densityOnDensityMatrix_, *X2, ONE);
    } else {
      Y2->elementWiseMultiply(ONE, *densityOnDensityMatrix_, *X2, ZERO);
    }
  }
  else
  {
    // First block row
    poissonOnPoissonMatrixOp_->apply(*X0, *Y0);
    poissonOnDensityMatrixOp_->apply(*X1, *Y0tmp);
    Y0->update(ONE, *Y0tmp, ONE);

    // Second block row
    if (hasDensityOnCms) {
      densityOnCmsMatrixOp_->apply(*X2, *Y1);
      Y1->elementWiseMultiply(ONE, *densityOnDensityMatrix_, *X1, ONE);
    } else {
      Y1->elementWiseMultiply(ONE, *densityOnDensityMatrix_, *X1, ZERO);
    }

    // Third block row
    cmsOnPoissonMatrixOp_->apply(*X0, *Y2);
    cmsOnDensityMatrixOp_->apply(*X1, *Y2tmp1);
    cmsOnCmsMatrixOp_->apply(*X2, *Y2tmp2);
    Y2->update(ONE, *Y2tmp1, ONE);
    Y2->update(ONE, *Y2tmp2, ONE);
  }
} //end Apply
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
Check
(bool verbose) const
{
  RCP<VEC > x = rcp(new VEC(getDomainMap()));
  RCP<VEC > b = rcp(new VEC(getRangeMap()));
  RCP<VEC > bb = rcp(new VEC(*b));
  x->randomize(); // Fill x with random numbers

  RCP<MV> x0 = x->offsetViewNonConst(poissonMap_, 0);
  RCP<MV> b0 = b->offsetViewNonConst(poissonMap_, 0);
  RCP<MV> x1 = x->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements());
  // Start x1 to view middle numDensity elements of x
  RCP<MV> b1 = b->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements());
  // Start b1 to view middle numDensity elements of b
  RCP<MV> b2 = b->offsetViewNonConst(cmsMap_, poissonMap_->getNodeNumElements()+densityMap_->getNodeNumElements());
  // Start b2 to view last numCms elements of b
  RCP<MV> x2 = x->offsetViewNonConst(cmsMap_, poissonMap_->getNodeNumElements()+densityMap_->getNodeNumElements());
  //Start x2 to view last numCms elements of x

  // The poisson-on-poisson matrix is singular. Make x0 orthogonal to the kernel.
  RCP<MV> ones = rcp(new MV(*b0));

  Scalar alpha;
  Array<Scalar> norm(1);
  Array<Scalar> dot(1);
  ArrayView<Scalar> norms = norm.view( 0, 1 );
  ArrayView<Scalar> dots = dot.view( 0, 1 );
  ones->putScalar( STS::one() );
  ones->norm2( norms );
  alpha = norms[0];
  alpha = STS::one() / alpha;
  ones->scale( alpha );
  ones->dot( *x0, dots );
  alpha = dots[0];
  alpha = -STS::one()*alpha;
  x0->update( alpha, *ones, STS::one() );

  apply(*x, *b); // Forward operation

  // Inverse if not exact, so we must modify b first:
  if (F_location_ == 1)
  {
    RCP<MV > DCx1 = rcp(new MV(*b2));
    densityOnCmsMatrixOp_->apply(*x1, *DCx1);
    b2->update(-STS::one(), *DCx1, STS::one()); // b2 = b2 - DC*x1
  }
  else
  {
    RCP<MV > DCx2 = rcp(new MV(*b1));
    densityOnCmsMatrixOp_->apply(*x2, *DCx2);
    b1->update(-STS::one(), *DCx2, STS::one()); // b1 = b1 - DC*x2
  }

  applyInverse(*b, *bb); // Reverse operation
  bb->update(-STS::one(), *x, STS::one()); // Should be zero
  Scalar resid = bb->norm2();

  if (verbose)
  {
    std::cout << "A22 self-check residual = " << resid << std::endl;
  } //end if

  TEUCHOS_TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");
} //end Check

TRAMONTO_INST_HELPER(dft_PolyA22_Coulomb_Tpetra_Operator)


//@HEADER
// ******************************************************************** 
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER

#ifndef DFT_POISSON_EPETRA_OPERATOR_H
#define DFT_POISSON_EPETRA_OPERATOR_H

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Ifpack.h"
#include "Teuchos_RefCountPtr.hpp"
#include <map>
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"

//! dft_Poisson_Epetra_Operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/    

class dft_Poisson_Epetra_Operator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors.
  //! Builds a Poisson operator for electrostatic systems.
  /*! Accepts matrix coefficients from the electrostatic operator of a Tramonto problem.  Can then later
      be used to apply that operator or apply its preconditioner.
      \param poissonMap (In) The global IDs of those equations that are part of the Poisson operator.
  */

  dft_Poisson_Epetra_Operator(const Epetra_Map & poissonMap);
  //@}
  //@{ \name Assembly methods.
  int initializeProblemValues();
  int insertMatrixValue(int rowGID, int colGID, double value);
  int finalizeProblemValues();
  //@}
  //@{ \name Destructor.
    //! Destructor
  ~dft_Poisson_Epetra_Operator();
  //@}
  
  //@{ \name Atribute set methods.

    //! Unsupported feature, returns -1.
  int SetUseTranspose(bool UseTranspose){return(-1);};
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a dft_Poisson_Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param X - (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y - (Out) An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of an inverse dft_Poisson_Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X - (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y - (Out) An Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
  */ 
  double NormInf() const {return(0.0);};
  //@}
  
  //@{ \name Atribute access functions

  //! Returns a character string describing the operator
  const char * Label() const{return(Label_);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(poissonMap_.Comm());};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(poissonMap_);};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(poissonMap_);};

  //! Returns a pointer to the Epetra_CrsMatrix object that is the Poisson matrix
  Epetra_CrsMatrix * getMatrix() { return(matrix_.get());}
  //@}
  
private:

  int insertRow();
  Epetra_Map poissonMap_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> matrix_;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> ifpackInverse_;
  char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool firstTime_;
  int curRow_;
  std::map<int, double> curRowValues_;
  Epetra_IntSerialDenseVector indices_;
  Epetra_SerialDenseVector values_;
  Ifpack factory_;
};

#endif /* DFT_POISSON_EPETRA_OPERATOR_H */
#include <iostream>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TMatrixDEigen.h>

using namespace std;

double det3By3(const TMatrixD A)
{
  assert( A.GetNrows()==3 && A.GetNcols()==3);
  double det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
    - A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] - A[2][2]*A[1][0]*A[0][1] ;
  return det;
}

TMatrixD rotationMatrix(int axis, double angle)
{
  TMatrixD r(3,3);
  r.Zero();
  if (axis!=0 && axis!=1 && axis!=2) return r;

  for( int i=0; i<3; i++ )
    {
      r[i][i]=cos(angle);
    }

  for( int i=-1; i<=1; i++ )
    {
      double row=(axis-i)%3; if(row<0) row+=3;
      double col=(axis+i)%3; if(col<0) col+=3;
      r[row][col]=i*sin(angle)+(1-i*i);
    }

  return r;
}

double dotProduct(const TVectorD a, const TVectorD b)
{
  assert( a.GetNrows() == 3 && b.GetNrows() == 3 );
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double dotProduct(const TMatrixD m, const TVectorD a)
{
  assert( m.GetNrows() == 3 && m.GetNcols() == 3 && a.GetNrows() == 3 );
  double ma0 = m[0][0]*a[0] + m[0][1]*a[1] + m[0][2]*a[2];
  double ma1 = m[1][0]*a[0] + m[1][1]*a[1] + m[1][2]*a[2];
  double ma2 = m[2][0]*a[0] + m[2][1]*a[1] + m[2][2]*a[2];
  double maArray[3] = {ma0,ma1,ma2};
  TVectorD ma(3,maArray);
  return dotProduct(a,ma);
}

double cofactor(const TMatrixD A, const int i, const int j)
{
  //cout << "Calculating cofactor " << i << "," << j << endl;
  if( A.GetNcols() != 3 || A.GetNrows() != 3 )
    {
      cout << "Matrix has the wrong dimensions!"  << endl;
      return 0.;
    }
  if( i < 0 || i >= 3 || j < 0 || j >= 3 )
    {
      cout << "Index out of bounds: " << i << " " << j << endl;
      return 0.;
    }

  //cout << "Input matrix:\n "
  //A.Print();
    
  int i0 = (i!=0) ? 0 : 1;
  int i1 = (i==2) ? 1 : 2;
  int j0 = (j!=0) ? 0 : 1;
  int j1 = (j==2) ? 1 : 2;

  int sign = ( (i+j)%2 == 0 ) ? 1 : -1;

  return sign*( A[i0][j0]*A[i1][j1] - A[i0][j1]*A[i1][j0] );
}

pair<TVectorD, TVectorD> factorDegenerateConic(const TMatrixD G, double zero=0.)
{
  //cout << "Calculating the decomposition of the ellipse into two lines" << endl;

  TVectorD Lplus(3), Lminus(3);
  Lplus.Zero();
  Lminus.Zero();

  if( G.GetNcols() != 3 || G.GetNrows() != 3 )
    {
      cout << "Matrix has the wrong dimensions!"  << endl;
      return make_pair(Lplus,Lminus);
    }

  if( G[0][0] == 0 && G[1][1] == 0 ) 
    {
      //corresponds to "horizontal and vertical" case in Table 1
      //cout << "Case: horizontal and vertical" << endl;
      Lplus.Zero();
      Lplus[0] = G[0][1];
      Lplus[2] = G[1][2];
      Lminus.Zero();
      Lminus[1] = G[0][1];
      Lminus[2] = G[0][2]-G[1][2];
      return make_pair(Lplus,Lminus);
    }

  TMatrixD Q(3,3);
  Q.Zero();

  bool swapXY = ( abs(G[0][0]) > abs(G[1][1]) ); //swap the x and y axes, ie, rotation around z axis by pi/2, then invert x axis

  if( swapXY )
    {
      //cout << "Swapping x and y axes" << endl;
      Q[0][0] = G[1][1]; Q[0][1] = G[1][0]; Q[0][2] = G[1][2];
      Q[1][0] = G[0][1]; Q[1][1] = G[0][0]; Q[1][2] = G[0][2];
      Q[2][0] = G[2][1]; Q[2][1] = G[2][0]; Q[2][2] = G[2][2];
      //check that this gives the same results
      //double piOver2 = 0.5*3.14159265359;
      //TMatrixD r1 = rotationMatrix(2,piOver2);
      //TMatrixD r2(3,3);
      //r2.Zero();
      //r2[0][0] = -1.;
      //r2[1][1] = 1.;
      //r2[2][2] = 1.;
      //Q = TMatrixD(r2,TMatrixD::kMult,TMatrixD(r1,TMatrixD::kMult,G));
    }
  else
    {
      Q = G;
    }

  Q *= (1./double(Q[1][1]));
  double q22 = cofactor(Q,2,2);

  if( q22 >= zero )
    {
      //corresponds to "parallel" case in Table 1
      //cout << "Case: parallel" << endl;
      double q00 = cofactor(Q,0,0);
      if( q00 > 0. )  
	{
	  //cout << "-q00 is negative!" << endl;
	  Lplus.Zero();
	  Lminus.Zero();
	  return make_pair(Lplus,Lminus);
	}
      Lplus.Zero();
      Lplus[0] = Q[0][1];
      Lplus[1] = Q[1][1];
      Lplus[2] = Q[1][2] - sqrt((-1.)*q00);
      Lminus.Zero();
      Lminus[0] = Q[0][1];
      Lminus[1] = Q[1][1];
      Lminus[2] = Q[1][2] + sqrt((-1.)*q00);
    }
  else
    {
      //corresponds to "intersecting" case in Table 1
      //cout << "Case: intersecting" << endl;
      double q02 = cofactor(Q,0,2);
      double q12 = cofactor(Q,1,2);
      if( q22 > 0. )
        {
          //cout << "-q22 is negative!" << endl;
          Lplus.Zero();
          Lminus.Zero();
          return make_pair(Lplus,Lminus);
        }
      double m = Q[0][1] - sqrt((-1.)*q22);
      Lplus.Zero();
      Lplus[0] = m;
      Lplus[1] = Q[1][1];
      Lplus[2] = (-1.)*Q[1][1]*(q12/q22) - m*(q02/q22);
      m = Q[0][1] + sqrt((-1.)*q22);
      Lminus.Zero();
      Lminus[0] = m;
      Lminus[1] = Q[1][1];
      Lminus[2] = (-1.)*Q[1][1]*(q12/q22) - m*(q02/q22);
    }

  if(swapXY)
    {
      // x <-> y
      double temp = Lplus[0];
      Lplus[0] = Lplus[1];
      Lplus[1] = temp;
      temp = Lminus[0];
      Lminus[0] = Lminus[1];
      Lminus[1] = temp;
    }

  return make_pair(Lplus,Lminus);
}

bool intersectionsOfEllipseAndLine(const TMatrixD ellipse, const TVectorD line, vector<TVectorD>& sols, double acceptLevel=1.e-12)
{
  //cout << "Calculating points of intersection between ellipse and line"

  sols.clear();

  if( ellipse.GetNcols() != 3 || ellipse.GetNrows() != 3 || line.GetNrows() !=3 )
    {
      cout << "Invalid size" << endl;
      return false;
    }

  double zeroArray[3] = {0.,0.,0.};
  TVectorD Zero(3,zeroArray);

  if( line == Zero )
    {
      //cout << "Empty line!" << endl;
      return false;
    }

  //Points of intersection between a line and a conic are eigenvectors of the cross product

  TMatrixD cross(3,3);
  cross[0][0] = line[1]*ellipse[2][0] - line[2]*ellipse[1][0];
  cross[0][1] = line[1]*ellipse[2][1] - line[2]*ellipse[1][1];
  cross[0][2] = line[1]*ellipse[2][2] - line[2]*ellipse[1][2];
  cross[1][0] = line[2]*ellipse[0][0] - line[0]*ellipse[2][0];
  cross[1][1] = line[2]*ellipse[0][1] - line[0]*ellipse[2][1];
  cross[1][2] = line[2]*ellipse[0][2] - line[0]*ellipse[2][2];
  cross[2][0] = line[0]*ellipse[1][0] - line[1]*ellipse[0][0];
  cross[2][1] = line[0]*ellipse[1][1] - line[1]*ellipse[0][1];
  cross[2][2] = line[0]*ellipse[1][2] - line[1]*ellipse[0][2];
  TMatrixDEigen intersect(cross);

  TMatrixD eigenvect = intersect.GetEigenVectors();

  //the 3d component is 1 in the homogeneous representation
  // => scale the matrix of eigenvectors
  for( int iCol=0; iCol<3; iCol++ )
    {
      for( int iRow=0; iRow<3; iRow++ )
	{
	  eigenvect[iRow][iCol] *= (1./double(eigenvect[2][iCol]));
	}
    }

  //if v is an eigenvector of LxM, then L.v=0 and vT.M.v=0
  //keep only the eigenvectors below the input accept level
  for( int iSol=0; iSol<3; iSol++ )
    {
      TVectorD thisSol(3);
      thisSol[0] = eigenvect[0][iSol];
      thisSol[1] = eigenvect[1][iSol];
      thisSol[2] = eigenvect[2][iSol];
      double dotProduct2 = pow(dotProduct(line,thisSol),2) + pow(dotProduct(ellipse,thisSol),2) ;
      if( dotProduct2 < acceptLevel ) 
	{
	  sols.push_back(thisSol);
	}
    }

  return true;
}

bool intersectionsOfEllipses(const TMatrixD A, const TMatrixD B, vector<TVectorD>& sols, double acceptLevel=1.e-6)
{
  //cout << "Calculating the points of intersection between two ellipses" << endl;

  //check matrix dimensions -- have to be 3x3
  if( A.GetNcols() != 3 || A.GetNrows() != 3 || B.GetNcols() != 3 || B.GetNrows() != 3 )
    {
      cout << "Matrix has the wrong dimensions!"  << endl;
      return false;
    }


  sols.clear();

  //the method will only work if A and B are invertible
  double detA = det3By3(A);
  double detB = det3By3(B);

  //cout << "The input ellipse matrices have determinant " << detA << " and " << detB << endl;

  if( detA == 0. || detB == 0. )
    {
      //cout << "One of the input matrices is not invertible!"  << endl;
      return false;
    }


  TMatrixD maxDet(3,3); 
  maxDet = ( abs(detB) > abs(detA) ) ? B : A;
  TMatrixD maxDetInv(maxDet);
  maxDetInv.Invert();

  TMatrixD minDet(3,3);
  minDet = ( abs(detB) > abs(detA) ) ? A : B;

  TMatrixD eigSyst(TMatrixD(maxDetInv,TMatrixD::kMult,minDet));
  TMatrixDEigen eig(eigSyst);

  TVectorD eigValRe = eig.GetEigenValuesRe();
  TVectorD eigValIm = eig.GetEigenValuesIm();
  int nRealEigenVal = 0;
  if( eigValRe.GetNrows() != 3 || eigValIm.GetNrows() != 3 )
    {
      cout << "Error calculating the eigenvalues!" << endl;
      return false;
    }

  int nIntersections = 0;
  double e;
  vector<TVectorD> intersections;

  for( int iEigenVal=0; iEigenVal<3; iEigenVal++ )
    { 
      if( abs(eigValIm[iEigenVal]) < acceptLevel ) 
	{
	  nRealEigenVal++;
	  e = eigValRe[iEigenVal];
	  pair<TVectorD, TVectorD> lines = factorDegenerateConic(TMatrixD(minDet-e*maxDet));
	  TVectorD Lplus  = lines.first;
	  TVectorD Lminus = lines.second;
	  //cout << "Checking for intersections of L+ and ellipse" << endl;
	  bool hasInt = intersectionsOfEllipseAndLine(maxDet,Lplus,intersections);
	  if( hasInt ) 
	    {
	      nIntersections += intersections.size();
	      for( int iSol=0 ; iSol<(int)intersections.size() ; iSol++ ) sols.push_back(intersections.at(iSol));
	    }
	  //cout << "Checking for intersections of L- and ellipse" << endl;
	  hasInt = intersectionsOfEllipseAndLine(maxDet,Lminus,intersections);
	  if( hasInt ) 
            {
              nIntersections += intersections.size();
	      for( int iSol=0 ; iSol<(int)intersections.size() ; iSol++ ) sols.push_back(intersections.at(iSol));
            }
	}
    }

  //cout << "There are " << nIntersections << " intersections between ellipses" << endl;
  //cout << "Intersections of the two ellipse:" << endl;
  //cout << "Neutrino momenta (Rochester method)" << endl;
  //for( int iSol=0 ; iSol<(int)nIntersections; iSol++ )
  //  {
  //    cout << "( " << sols.at(iSol)[0] << " , " << sols.at(iSol)[1] << " , " << sols.at(iSol)[2] << " )" << endl;
  //  }
  
  if(nIntersections==0) return false;
  return true;
}

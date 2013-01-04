#include "SpectralCluster.h"

namespace statistics
{

void SpectralCluster::ConvertToSimpleEigenvalue()
{
  std::vector<double> degree;
  degree.resize(similarity_matrix_.size(), 0);
  for(unsigned int i = 0; i < degree.size(); ++i)
    for(unsigned int j = 0; j < similarity_matrix_[i].size(); ++j)
      degree[i]+= similarity_matrix_[i][j];

  for(unsigned int i = 0; i < degree.size(); ++i)
    for(unsigned int j = 0; j < similarity_matrix_[i].size(); ++j)
    {
      if( i == j)
        similarity_matrix_[i][j] = degree[i] - similarity_matrix_[i][j];
      else
        similarity_matrix_[i][j] = - similarity_matrix_[i][j];

      similarity_matrix_[i][j] = similarity_matrix_[i][j] * 
          (1 / std::sqrt(degree[i])) * (1 / std::sqrt(degree[j]));
    }

}

bool SpectralCluster::PerformEigenAnalysis(
    std::vector<std::vector<double> > &vectors, std::vector<double> &values)
{
  utilities::CreateIdentityMatrix(vectors, similarity_matrix_.size(), 
      similarity_matrix_.size());
  values.resize( similarity_matrix_.size() );
  Schur(vectors, values);
  SortEigenvalues(vectors, values);
  return true;
}

bool SpectralCluster::FindNextPivot(std::vector<std::vector<double> > &M, 
    unsigned int &row, unsigned int &col)
{
  double maxval = 0;
  bool found = false;
  for(unsigned int r = 0; r < M.size(); ++r)
    for(unsigned int c = 0; c < M[r].size(); ++c)
      //if(r != c && std::fabs(M[r][c]) > (EPS_ * std::fabs(M[r][r] - M[c][c])))
      if(r != c && std::fabs(M[r][c]) > EPS_)
      {
        if(std::fabs(M[r][r] - M[c][c]) < EPS_) // Nearly identical
        {
          row = r;
          col = c;
          return true;
        }
        else if( (std::fabs(M[r][c]) / std::fabs(M[r][r] - M[c][c])) > maxval )
        {
          row = r;
          col = c;
          maxval = (std::fabs(M[r][c]) / std::fabs(M[r][r] - M[c][c]));
          found = true;
        }
      }
  return found;
}

bool SpectralCluster::Schur(std::vector<std::vector<double> > &U, 
    std::vector<double> &e)
{
  std::vector<std::vector<double> > M = similarity_matrix_;
  std::vector<std::vector<double> > Usub;
  int iter = 0;
  int itermax = 500;
  unsigned int row, col;

  utilities::CreateIdentityMatrix(Usub, 2, 2);
  while(FindNextPivot(M, row, col) && (iter < itermax) )
  {
    iter++;
    std::cout<<iter<<" "<<row<<" "<<col<<" "<<M[row][col]<<std::endl;
    double a, b, c, L;
    //if( std::fabs(M[row][col]) > (std::fabs(M[row][row] - M[col][col])) )
    if(EPS_ < std::fabs(M[row][row] - M[col][col]))
    {
      a = 2 * M[row][col] / (M[row][row] - M[col][col]);
      b = 2 * M[col][row] / (M[row][row] - M[col][col]);
      c = std::sqrt(1 + (a*b));
      L = std::sqrt(std::pow(a, 2) + std::pow(1+c, 2));
      Usub[0][0] = (1 + c) / L; Usub[0][1] = -a / L;
      Usub[1][0] = a / L;       Usub[1][1] = (1+c) / L;
    }
    else
    {
      a = std::fabs(M[row][col]);
      b = std::fabs(M[col][row]);
      L = std::sqrt(a + b);
      Usub[0][0] = std::sqrt(b) / L; Usub[0][1] = -std::sqrt(a) / L;
      Usub[1][0] = std::sqrt(a) / L; Usub[1][1] = std::sqrt(b) / L;
    }
    for(unsigned int r = 0; r < M.size(); ++r)
    {
      double newr = M[r][row]*Usub[0][0] + M[r][col]*Usub[1][0];
      double newc = M[r][row]*Usub[0][1] + M[r][col]*Usub[1][1];
      M[r][row] = newr;
      M[r][col] = newc;
    }
    for(unsigned int c = 0; c < M[0].size(); ++c)
    {
      double newr = M[row][c]*Usub[0][0] + M[col][c]*Usub[1][0];
      double newc = M[row][c]*Usub[0][1] + M[col][c]*Usub[1][1];
      M[row][c] = newr;
      M[col][c] = newc;
    }
    for(unsigned int r = 0; r < M.size(); ++r)
    {
      double newr = U[r][row]*Usub[0][0] + U[r][col]*Usub[1][0];
      double newc = U[r][row]*Usub[0][1] + U[r][col]*Usub[1][1];
      U[r][row] = newr;
      U[r][col] = newc;
    }
  }

  e = utilities::MatrixDiagonal(M);
  return true;
}

bool SpectralCluster::Jacobi2(std::vector<std::vector<double> > &U, 
    std::vector<double> &e)
{
  std::vector<std::vector<double> > A = similarity_matrix_;
  unsigned int size = A.size();
  e = utilities::MatrixDiagonal(A);
  std::vector<double> bw = e;
  int iter = 0;
  int itermax = 10000;

  for(unsigned int i = 0; i < size; ++i)
    A[i][i] = 0;

  while( iter < itermax)
  {
    double tresh=utilities::FrobeniusNorm(utilities::UpperTriangularMatrix(A)) /
        (size*4);
    if( tresh < EPS_)
      return true;
    std::vector<double> zw(size, 0);
    for(unsigned int p = 0; p < size; ++p)
    {
      for(unsigned int q = 0; q < size; ++q)
      {
        double gapq = 10 * std::fabs(A[p][q]);
        double termp = gapq + std::fabs( e[p] );
        double termq = gapq + std::fabs( e[q] );
        double t = 0, c = 0, s = 0, h = 0, g = 0;
        double theta = 0, tau = 0;
        if( iter > 3 && termp == std::fabs(e[p]) && termq == std::fabs(e[q]))
        {
          A[p][q] = 0; // Annihilate tiny elements
        }
        else
        {
          if( std::fabs(A[p][q]) >= tresh)
          {
            h = e[q] - e[p];
            double term= std::fabs(h) + gapq;
            if( term == std::fabs(h) )
            {
              t = A[p][q] / h;
            }
            else
            {
              theta = (0.5*h) / A[p][q];
              t = 1 / (std::fabs(theta) + std::sqrt(1 + std::pow(theta,2)));
              if( theta < 0)
                t = -t;
            }
            c = 1 / std::sqrt(1 + std::pow(t, 2));
            s = t * c;
            tau = s / (1 + c);
            h = t * A[p][q];
            zw[p] = zw[p] - h;
            zw[q] = zw[q] + h;
            e[p] = e[p] - h;
            e[q] = e[q] + h;
            A[p][q] = 0;
            for(unsigned int j = 0; j < p; ++j)
            {
              g = A[j][p];
              h = A[j][q];
              A[j][p] = g - (s * (h + (g*tau)));
              A[j][q] = h + (s * (g - (h*tau)));
            }
            for(unsigned int j = p+1; j < q; ++j)
            {
              g = A[p][j];
              h = A[j][q];
              A[p][j] = g - (s * (h + (g*tau)));
              A[j][q] = h + (s * (g - (h*tau)));
            }
            for(unsigned int j = q+1; j < size; ++j)
            {
              g = A[p][j];
              h = A[q][j];
              A[p][j] = g - (s * (h + (g*tau)));
              A[q][j] = h + (s * (g - (h*tau)));
            }
            for(unsigned int j = 0; j < size; ++j)
            {
              g = U[j][p];
              h = U[j][q];
              U[j][p] = g - (s * (h + (g*tau)));
              U[j][q] = h + (s * (g - (h*tau)));
            }
          }
        }
      }
    }
    for(unsigned int i = 0; i < bw.size(); ++i)
      bw[i] += zw[i];
    e=bw;
    ++iter;
  }
  return true;
}

bool SpectralCluster::Jacobi(std::vector<std::vector<double> > &vectors, 
    std::vector<double> &values)
{
  unsigned int col, row;
  bool iterating = true;
  unsigned int n = similarity_matrix_.size();
  double previous_max = 0;

  std::vector<std::vector<double> > M, S;
  S = similarity_matrix_;
 
  while(iterating)
  {
    M = S;
    utilities::MatrixAbs(M);

    for (unsigned int k=0; k<n; k++)
      M[k][k]=0;

    utilities::MaxElementInMatrix(M, row, col);
    double Smax = S[row][col];
    std::cout<<Smax<<" "<<row<<" "<<col<<std::endl;
    if (row == col)
    {
      for (unsigned int i=0; i<n; i++) 
        values[i] = S[i][i];
      return true;
    }
    ZeroRowCol(S, vectors, row, col);
    if (std::fabs(Smax - previous_max) < EPS_ )//* utilities::FrobeniusNorm(S)) 
      iterating = false; 

    previous_max = Smax;
  }
  for (unsigned int i=0; i<n; i++) 
    values[i] = S[i][i];

  return true;
}

void SpectralCluster::ZeroRowCol(std::vector<std::vector<double> > &S, 
    std::vector<std::vector<double> > &vectors, unsigned int row, 
    unsigned int col)
{
  double t, c, s, theta;
  unsigned int n = S.size();

  if (row == col) 
    return;

  theta=(S[row][row]-S[col][col])/(2*S[row][col]);
  if (theta < EPS_)
    t = 1/(abs(theta)+sqrt(theta*theta+1));
  else 
    t = 1/abs(2*theta);

  if (theta<0) 
    t = -t;

  c = 1/std::sqrt(t*t + 1);
  s = c*t;
  std::vector<std::vector<double> >R, Rprime;
  utilities::CreateIdentityMatrix(R, n, n);
  R[row][row] = c;
  R[col][col] = c; 
  R[row][col] = s; 
  R[col][row] = -s;
  Rprime = utilities::TransposeMatrix(R);
  S = utilities::MatrixProduct(S, Rprime);
  S = utilities::MatrixProduct(R, S);
  vectors = utilities::MatrixProduct(vectors, R); 
  return;
}

void SpectralCluster::SortEigenvalues(std::vector<std::vector<double> >&vectors,
    std::vector<double> &values)
{
  std::vector< std::pair<double, int> > e;
  std::vector<std::vector<double> > U = vectors;
  e.resize( values.size() );
  for(unsigned int i = 0; i < values.size(); ++i)
  {
    e[i].first = values[i];
    e[i].second = i;
  }
  std::sort(e.begin(), e.end());
  std::reverse(e.begin(), e.end());
  for(unsigned int i = 0; i < values.size(); ++i)
  {
    values[i] = e[i].first;
    for(unsigned int j = 0; j < values.size(); ++j)
      U[j][i] = vectors[j][e[i].second];
  }
  vectors = U;
}

}

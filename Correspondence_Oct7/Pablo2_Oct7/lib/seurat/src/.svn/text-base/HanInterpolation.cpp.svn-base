#include "HanInterpolation.h"



HanInterpolation::HanInterpolation(){
	m_M3DQuadFigure = 0;
	m_InterpolationLevel = 0;
}
HanInterpolation::~HanInterpolation(){
	ClearVectors();
}


void HanInterpolation::Update(){

	if(m_M3DQuadFigure==0){
#if WIN32
		throw exception("QuadFigure not set in HanInterpolation");;
#else
        throw exception();
#endif
	}

	ClearVectors();	
	CalculateSradSeAndDerivatives();
	InterpolatePoints();	
}


void HanInterpolation::CalculateSradSeAndDerivatives(){
//Calculation of rSrad and rSe for all the spokes
	for(int i = 0; i < m_M3DQuadFigure->getRowCount(); i++){
		
		vector< vector< Matrix > > vvMatrix;
		m_rSradMatrix.push_back(vvMatrix);

		vector< vector< Vector3D > > vvDerivatives;
		m_Derivatives.push_back(vvDerivatives);

		for(int j = 0; j < m_M3DQuadFigure->getColumnCount(); j++){	

			vector< Matrix > vMatrix;
			m_rSradMatrix[i].push_back(vMatrix);

			vector< Vector3D > vDerivatives;
			m_Derivatives[i].push_back(vDerivatives);

			//     prim2   _____    prim3
			//		  u2 /	    /
			//			/__u1_ /
			//  prim0		     prim1
			M3DPrimitive* prim0 = m_M3DQuadFigure->getPrimitivePtr(i,j);
			M3DPrimitive* prim1 = 0;
			M3DPrimitive* prim2 = 0;

			Vector3D du1;
			Vector3D du2;
			Vector3D dUdu1;
			Vector3D dUdu2;
			Vector3D U0;
			double r0 = 0;

			//compute the derivatives in u1 and u2 direction at current position
			GetPDerivatives(i,j,du1,du2, &prim1, &prim2);
						

			m_Derivatives[i][j].push_back(du1);
			m_Derivatives[i][j].push_back(du2);

			

			//compute the derivatives of the U vectors according to the side
			for(int side = 0; side < 2; side++){
				
				if(side == 0){
					r0 = dynamic_cast<M3DQuadPrimitive*>(prim0)->getR0();
					
					U0 = prim0->getU0();

					if(i < m_M3DQuadFigure->getRowCount() - 1){
						if(i == 0){							
							dUdu1 = (prim1->getU0()) - (U0);
						}else{
							M3DPrimitive* tempprim0 = m_M3DQuadFigure->getPrimitivePtr(i-1,j);
							dUdu1 = ((prim1->getU0()) - (tempprim0->getU0()))/2.0;
						}
					}else{
						dUdu1 = (prim0->getU0()) - (prim1->getU0());
					}
					
					if(j < m_M3DQuadFigure->getColumnCount() - 1){
						if(j == 0){
							dUdu2 = (prim2->getU0()) - (U0);
						}else{
							M3DPrimitive* tempprim0 = m_M3DQuadFigure->getPrimitivePtr(i,j-1);
							double r0 = dynamic_cast<M3DQuadPrimitive*>(tempprim0)->getR0();
							dUdu2 = (prim2->getU0() - tempprim0->getU0())/2.0;
						}
					}else{
						dUdu2 = U0 - prim2->getU0();
					}										
				}else{
					double r0 = dynamic_cast<M3DQuadPrimitive*>(prim0)->getR1();
					
					U0 = prim0->getU1();

					if(i+1 < m_M3DQuadFigure->getRowCount()){
						if(i == 0){							
							dUdu1 = prim1->getU1() - U0;
						}else{
							M3DPrimitive* tempprim0 = m_M3DQuadFigure->getPrimitivePtr(i-1,j);							
							dUdu1 = (prim1->getU1() - tempprim0->getU1())/2.0;
						}
					}else{
						dUdu1 = U0 - prim1->getU1();
					}
					
					if(j < m_M3DQuadFigure->getColumnCount() - 1){
						if(j == 0){
							dUdu2 = prim2->getU1() - U0;
						}else{
							M3DPrimitive* tempprim0 = m_M3DQuadFigure->getPrimitivePtr(i,j-1);
							dUdu2 = (prim2->getU1() - tempprim0->getU1())/2.0;
						}
					}else{
						dUdu2 = U0 - prim2->getU1();
					}										
				}



				//store the rSrad matrices calculated according to the side
				Matrix rSrad;
				GetRSRadMatrix(du1, du2, dUdu1, dUdu2, U0, r0, rSrad);
				m_rSradMatrix[i][j].push_back(rSrad);					
			}


			if(i==0 || i == m_M3DQuadFigure->getRowCount() - 1
				|| j == 0 || j == m_M3DQuadFigure->getColumnCount() - 1){


				Vector3D du;
				Vector3D dndu;
				Vector3D dUdu;
				GetEdgeDerivatives(i,j,du,dndu,dUdu);

				double r0 = dynamic_cast<M3DQuadEndPrimitive*>(prim0)->getREnd();					
				Vector3D U0 = dynamic_cast<M3DQuadEndPrimitive*>(prim0)->getUEnd();

				Matrix rSe;
				GetRSEMatrix(du, dndu, dUdu, U0, r0, rSe);
				m_rSradMatrix[i][j].push_back(rSe);

			}


		}
	}
}

void HanInterpolation::InterpolatePoints(){

	int rowsize =  m_M3DQuadFigure->getRowCount() + (pow(2.0, m_InterpolationLevel) - 1)*( m_M3DQuadFigure->getRowCount()-1);
	//int colsize =  m_M3DQuadFigure->getColumnCount() + (pow(2.0, m_InterpolationLevel) - 1)*( m_M3DQuadFigure->getColumnCount()-1);
	
	vector< Vector3D > mCol;
	vector< vector< Vector3D > > mRow;
	mRow.insert(mRow.begin(), rowsize, mCol);
		
	m_InterpolatedPoints.insert(m_InterpolatedPoints.begin(), 2, mRow);
	m_InterpolatedSpokes.insert(m_InterpolatedSpokes.begin(), 2, mRow);

	for(int i = 0; i < m_M3DQuadFigure->getRowCount() - 1; i++){
		for(int j = 0; j < m_M3DQuadFigure->getColumnCount() - 1; j++){
			

			//     prim2   _____    prim3
			//		  u2 /	    /
			//			/__u1_ /
			//  prim0		     prim1
			M3DPrimitive* prim0 = m_M3DQuadFigure->getPrimitivePtr(i,j);
			M3DPrimitive* prim1 = m_M3DQuadFigure->getPrimitivePtr(i+1,j);
			M3DPrimitive* prim2 = m_M3DQuadFigure->getPrimitivePtr(i,j+1);
			M3DPrimitive* prim3 = m_M3DQuadFigure->getPrimitivePtr(i+1,j+1);

			Vector3D p0 = prim0->getX();
			Vector3D p1 = prim1->getX();
			Vector3D p2 = prim2->getX();
			Vector3D p3 = prim3->getX();

			vector< Vector3D > dp0 = m_Derivatives[i][j];
			vector< Vector3D > dp1 = m_Derivatives[i+1][j];
			vector< Vector3D > dp2 = m_Derivatives[i][j+1];
			vector< Vector3D > dp3 = m_Derivatives[i+1][j+1];

			//normals to the medial sheet
			Vector3D n0 = prim0->getU1() - prim0->getU0();
			n0.normalize();
			Vector3D n1 = prim1->getU1() - prim1->getU0();
			n1.normalize();
			Vector3D n2 = prim2->getU1() - prim2->getU0();
			n2.normalize();
			Vector3D n3 = prim3->getU1() - prim3->getU0();
			n3.normalize();


			Vector3D Hm[4][4];//hermite matrix for the current quad on the medial sheet
			GetHermiteMatrix(p0, p1, p2, p3,
							 dp0,dp1,dp2,dp3,
							 n0, n1, n2, n3, Hm);

			
			//This interpolates the points in the quad
			InterpolateMedial(i, j, prim0, Hm);

			if(j == 0){				
				InterpolateCrest(i, 0, i + 1, 0, Hm);
			}else if(j == m_M3DQuadFigure->getColumnCount()-2){
				InterpolateCrest(i, m_M3DQuadFigure->getColumnCount()-1, i+1, m_M3DQuadFigure->getColumnCount()-1, Hm);
			}

			if(i==0){		
				InterpolateCrest(0, j, 0, j+1, Hm);
			}else if(i == m_M3DQuadFigure->getRowCount()-2){
				InterpolateCrest(m_M3DQuadFigure->getRowCount()-1, j, m_M3DQuadFigure->getRowCount()-1, j + 1, Hm);		
			}
		}
	}
}

void HanInterpolation::InterpolateCrest(int i,int j,int inext, int jnext, Vector3D Hm[4][4]){
	
	M3DPrimitive* prim0 = m_M3DQuadFigure->getPrimitivePtr(i,j);

	Vector3D U = dynamic_cast<M3DQuadEndPrimitive*>(prim0)->getUEnd();
	double r = dynamic_cast<M3DQuadEndPrimitive*>(prim0)->getREnd();
	Vector3D n = prim0->getU1() - prim0->getU0();
	

	Matrix rSe0 = m_rSradMatrix[i][j][2];
	Matrix rSe1 = m_rSradMatrix[inext][jnext][2];

	//compute eigenvalues and eigenvectors from sRad matrices
	Vector L1, L2;
	Matrix V1, V2;

	//cout<<rSrad0(0,0)<<", "<<rSrad0(0,1)<<";"<<rSrad0(1,0)<<","<<rSrad0(1,1)<<endl;
	rSe0.factorEV(L1,V1, NON_SYM);
	//cout<<"e= "<<L11(0)<<" "<<L11(1)<<endl;
	//cout<<V11(0,0)<<" "<<V11(0,1)<<endl;
	//cout<<V11(1,0)<<" "<<V11(1,1)<<endl;								

	rSe1.factorEV(L2,V2, NON_SYM);

	double lam11, lam12, lam21, lam22;
	Vector e11, e12, e21, e22;

	//ordering eigenvalues and eigenvectors
	int a = 0;
	int b = 1;
	if (L1(a) < L1(b)){
		a = 1;
		b = 0;
	}
	
	lam11 = L1(a);
	lam12 = L1(b);

	e11 = V1.getColumn(a);
	e12 = V1.getColumn(b);

	a = 0;
	b = 1;
	if (L2(a) < L2(b)){
		a = 1;
		b = 0;
	}
	
	lam21 = L2(a);
	lam22 = L2(b);

	e21 = V2.getColumn(a);
	e22 = V2.getColumn(b);	

	double logLam11 = GetLogLambda(lam11);
	double logLam12 = GetLogLambda(lam12);
	double logLam21 = GetLogLambda(lam21);
	double logLam22 = GetLogLambda(lam22);

	double step = pow((double)2, (double)m_InterpolationLevel);
	double stepsize = 1/step;

	//rowpos is the global position in the grid of the interpolated points
	for(double u = 0; u < 1; u+=stepsize){

		double stepu = u/10;
		Vector3D spokeintegral(0,0,0);
		Vector3D x0 = prim0->getX();
		
		if(stepu != 0){
			for(double curru = 0; curru <= u; curru+=stepu){
				//weighted average
				// To find the new rSe, need to calculate lambdas and eigenvectors
				double logAvg1 = (1-curru)*logLam11 + (curru)*logLam21;
				double logAvg2 = (1-curru)*logLam12 + (curru)*logLam22;

				double lam1 = 1 - exp(logAvg1);
				double lam2 = 1 - exp(logAvg2);

				double avgx1, avgy1, avgx2, avgy2, avgtheta1, avgtheta2;

				avgx1 = (1-curru)*e11(0) + (curru)*e21(0);
				avgy1 = (1-curru)*e11(1) + (curru)*e21(1);
				avgtheta1 = atan2(avgy1, avgx1);							

				avgx2 = (1-curru)*e12(0) + (curru)*e22(0);
				avgy2 = (1-curru)*e12(1) + (curru)*e22(1);				
				avgtheta2 = atan2(avgy2, avgx2);						

				double neweigenv[4] = { cos(avgtheta1), cos(avgtheta2), sin(avgtheta1), sin(avgtheta2) };

				Matrix NewV = Matrix(2,2,neweigenv,true);

				double newlambda[4] = { lam1, 0, 0, lam2 };

				Matrix NewL = Matrix(2,2,newlambda,true);		
				
				
				Matrix NewVi;
				NewV.inverse(NewVi);

				Matrix NewrSe = NewV * NewL * NewVi;

				Vector3D x1;				

				if(i - inext == 0){
					GetInterpolatedPoint(0,curru,Hm, x0);
					GetInterpolatedPoint(0,curru+stepu,Hm, x1);
				}else{
					GetInterpolatedPoint(curru,0,Hm, x0);
					GetInterpolatedPoint(curru+stepu,0,Hm, x1);
				}

				Matrix dEdun(2,3);
				GetdEdu(x0,x1,n,U, NewrSe, dEdun );

				Vector dEdu = dEdun.getRow(0);
				Vector dEdn = dEdun.getRow(1);

				Vector3D dEu = Vector3D(dEdu(0), dEdu(1), dEdu(2));
				Vector3D dEn = Vector3D(dEdn(0), dEdn(1), dEdn(2));
				
				spokeintegral += (stepu * dEu) + (stepu * dEn);
			}
		}

		Vector3D newSpoke = U*r + spokeintegral;

	}
}

void HanInterpolation::GetdEdu(Vector3D x0, Vector3D x1, Vector3D n, Vector3D U0, Matrix rSe, Matrix& dEdun){
	Vector3D pu1;
	pu1 = x1 - x0;

	Matrix pun(2, 3, pu1.getX(), pu1.getY(), pu1.getZ(), n.getX(), n.getY(), n.getZ());

	Matrix U(1,3, U0.getX(), U0.getY(), U0.getZ());

	//this is Damon(2004) compatibility condition
	//implies that a spoke S is perpedicular to the implied boundary (2.12)
	Matrix rUn = -1 * pun * U.t();
	
	Matrix I(3, 3, 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);//identity matrix
	
	Matrix Q = pun * ((U.t() * U) - I);
	
	dEdun = (rSe.t() * Q) + (rUn * U);
}

void HanInterpolation::InterpolateMedial(int i, int j, M3DPrimitive* prim0, Vector3D Hm[4][4]){
	for(int side = 0; side < 2; side++){
								
		//Jared's interpolation is done here

		Matrix rSrad0 = m_rSradMatrix[i][j][side];
		Matrix rSrad1 = m_rSradMatrix[i+1][j][side];
		Matrix rSrad2 = m_rSradMatrix[i][j+1][side];
		Matrix rSrad3 = m_rSradMatrix[i+1][j+1][side];				
	
		//compute eigenvalues and eigenvectors from sRad matrices
		Vector L11, L21, L12, L22;
		Matrix V11, V21, V12, V22;

		//cout<<rSrad0(0,0)<<", "<<rSrad0(0,1)<<";"<<rSrad0(1,0)<<","<<rSrad0(1,1)<<endl;
		rSrad0.factorEV(L11,V11, NON_SYM);
		//cout<<"e= "<<L11(0)<<" "<<L11(1)<<endl;
		//cout<<V11(0,0)<<" "<<V11(0,1)<<endl;
		//cout<<V11(1,0)<<" "<<V11(1,1)<<endl;								

		rSrad1.factorEV(L21,V21, NON_SYM);
		rSrad2.factorEV(L12,V12, NON_SYM);
		rSrad3.factorEV(L22,V22, NON_SYM);

		double Lam1_11, Lam2_11, Lam1_21, Lam2_21, Lam1_12, Lam2_12, Lam1_22, Lam2_22;
		Vector e1_11, e2_11, e1_21, e2_21, e1_12, e2_12, e1_22, e2_22;

		//ordering eigenvalues and eigenvectors
		if (L11(0) > L11(1))
		{
			Lam1_11 = L11(0);
			Lam2_11 = L11(1);

			e1_11 = V11.getColumn(0);
			e2_11 = V11.getColumn(1);
		}
		else
		{
			Lam1_11 = L11(1);
			Lam2_11 = L11(0);

			e1_11 = V11.getColumn(1);
			e2_11 = V11.getColumn(0);
		}

		// Corner 2
		if (L21(0) > L21(1))
		{
			Lam1_21 = L21(0);
			Lam2_21 = L21(1);

			e1_21 = V21.getColumn(0);
			e2_21 = V21.getColumn(1);
		}
		else
		{
			Lam1_21 = L21(1);
			Lam2_21 = L21(0);

			e1_21 = V21.getColumn(1);
			e2_21 = V21.getColumn(0);
		}

		//Corner 3
		if (L12(0) > L12(1))
		{
			Lam1_12 = L12(0);
			Lam2_12 = L12(1);

			e1_12 = V12.getColumn(0);
			e2_12 = V12.getColumn(1);
		}
		else
		{
			Lam1_12 = L12(1);
			Lam2_12 = L12(0);

			e1_12 = V12.getColumn(1);
			e2_12 = V12.getColumn(0);
		}

		//Corner 4
		if (L22(0) > L22(1))
		{
			Lam1_22 = L22(0);
			Lam2_22 = L22(1);

			e1_22 = V22.getColumn(0);
			e2_22 = V22.getColumn(1);
		}
		else
		{
			Lam1_22 = L22(1);
			Lam2_22 = L22(0);

			e1_22 = V22.getColumn(1);
			e2_22 = V22.getColumn(0);
		}

		double logLam1_11 = GetLogLambda(Lam1_11);
		double logLam2_11 = GetLogLambda(Lam2_11);
		double logLam1_21 = GetLogLambda(Lam1_21);
		double logLam2_21 = GetLogLambda(Lam2_21);
		double logLam1_12 = GetLogLambda(Lam1_12);
		double logLam2_12 = GetLogLambda(Lam2_12);
		double logLam1_22 = GetLogLambda(Lam1_22);
		double logLam2_22 = GetLogLambda(Lam2_22);				

		Vector3D U;
		double r=0;
		if(side == 0){
			U = prim0->getU0();
			M3DQuadPrimitive* tempprim = dynamic_cast<M3DQuadPrimitive*>(prim0);
			r = tempprim->getR0();
		}else{
			U = prim0->getU1();
			M3DQuadPrimitive* tempprim = dynamic_cast<M3DQuadPrimitive*>(prim0);
			r = tempprim->getR1();
		}

		
		Vector3D x0 = prim0->getX();//interpolated medial point, the initial point is at prim0
		Vector3D spokeintegral(0,0,0);//integration over a straight line to new interpolated spoke

		double step = pow((double)2, (double)m_InterpolationLevel);
		double stepsize = 1/step;

		double limu = 1;
		if(i == m_M3DQuadFigure->getRowCount() - 2){//If is the last row increment the limit to take care of it
			limu += stepsize;//this is done in order to iterpolate over the last row
		}
		double limv = 1;
		if(j == m_M3DQuadFigure->getColumnCount() - 2){//if its the last column increment the limit to take care of it
			limv += stepsize;//this is done in order to iterpolate over the last column
		}

		

		//rowpos is the global position in the grid of the interpolated points
		for(double u = 0, rowpos = i*step; u < limu; u+=stepsize, rowpos++){
			for(double v = 0; v < limv; v+=stepsize){	

				double stepu = u/10;
				double stepv = v/10;						
				
				if(stepu != 0 || stepv != 0){
					for(double curru = 0, currv = 0; curru <= u && currv <= v; curru+=stepu, currv+=stepv){
						//weighted average
						// To find rSrad, need to calculate lambdas and eigenvectors
						double logAvg1 = (1-curru)*(1-currv)*logLam1_11 + (curru)*(1-currv)*logLam1_21 + (1-curru)*(currv)*logLam1_12 + (curru)*(currv)*logLam1_22;
						double logAvg2 = (1-curru)*(1-currv)*logLam2_11 + (curru)*(1-currv)*logLam2_21 + (1-curru)*(currv)*logLam2_12 + (curru)*(currv)*logLam2_22;
					
						double Lam1 = 1 - exp(logAvg1);
						double Lam2 = 1 - exp(logAvg2);

						double avgx1, avgy1, avgx2, avgy2, avgtheta1, avgtheta2;

						avgx1 = (1-curru)*(1-currv)*e1_11(0) + (curru)*(1-currv)*e1_21(0) + (1-curru)*(currv)*e1_12(0) + (curru)*(currv)*e1_22(0);
						avgy1 = (1-curru)*(1-currv)*e1_11(1) + (curru)*(1-currv)*e1_21(1) + (1-curru)*(currv)*e1_12(1) + (curru)*(currv)*e1_22(1);
						avgtheta1 = atan2(avgy1, avgx1);							

						avgx2 = (1-u)*(1-v)*e2_11(0) + (u)*(1-v)*e2_21(0) + (1-u)*(v)*e2_12(0) + (u)*(v)*e2_22(0);
						avgy2 = (1-u)*(1-v)*e2_11(1) + (u)*(1-v)*e2_21(1) + (1-u)*(v)*e2_12(1) + (u)*(v)*e2_22(1);
						avgtheta2 = atan2(avgy2, avgx2);						

						//double neweigenv[4] = { cos(avgtheta1), sin(avgtheta1), cos(avgtheta2), sin(avgtheta2) }; The values must be transposed 
						//P = (er1(u1, u2)^T er2(u1, u2)^T
						double neweigenv[4] = { cos(avgtheta1), cos(avgtheta2),sin(avgtheta1), sin(avgtheta2) };

						Matrix NewV = Matrix(2,2,neweigenv,true);

						double newlambda[4] = { Lam1, 0, 0, Lam2 };

						Matrix NewL = Matrix(2,2,newlambda,true);		
						
						
						Matrix NewVi;
						NewV.inverse(NewVi);

						Matrix NewrSrad = NewV * NewL * NewVi;								
						
						GetInterpolatedPoint(curru,currv,Hm, x0);						

						Vector3D x1;
						if(curru < 1){
							GetInterpolatedPoint(curru+stepu,currv,Hm, x1);
						}else{
							GetInterpolatedPoint(curru-stepu,currv,Hm, x1);
						}
						Vector3D x2;
						if(v < 1){
							GetInterpolatedPoint(curru,currv+stepv,Hm, x2);
						}else{
							GetInterpolatedPoint(curru,currv-stepv,Hm, x2);
						}
						
						Matrix dSdun(2,3);
						GetdSdu(x0,x1,x2,stepu,stepv,U, NewrSrad, dSdun );

						Vector dSdu = dSdun.getRow(0);
						Vector dSdv = dSdun.getRow(1);

						Vector3D dSu = Vector3D(dSdu(0), dSdu(1), dSdu(2));
						Vector3D dSv = Vector3D(dSdv(0), dSdv(1), dSdv(2));
						
						spokeintegral += (stepu * dSu) + (stepv * dSv);

					}
				}

				Vector3D newSpoke = U*r + spokeintegral;
				
				m_InterpolatedPoints[side][rowpos].push_back(x0);
				m_InterpolatedSpokes[side][rowpos].push_back(newSpoke);
			}
		}							
	}
}

void HanInterpolation::ClearVectors(){

	for(int i = 0; i < m_rSradMatrix.size(); i++){
		for(int j = 0; j < m_rSradMatrix[i].size(); j++){
			m_rSradMatrix[i][j].clear();
			m_Derivatives[i][j].clear();
		}
		m_rSradMatrix[i].clear();
		m_Derivatives[i].clear();
	}
		
	m_rSradMatrix.clear();
	m_Derivatives.clear();
	m_InterpolatedPoints.clear();
	m_InterpolatedPoints.clear();

}


/*
*   \fn Vector3D* getHermiteMatrix(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3, Vector3D n0, Vector3D n1, Vector3D n2, Vector3D n3)
*	\brief Calculates the corresponding hermiteMatrix from the 4 control points and the 4 given normals
*	\param Vector3D p0-p1 the control points
*	\param Vector3D n0-n1 the normals of the control points
*/
void HanInterpolation::GetHermiteMatrix(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3, vector< Vector3D > dp0, vector< Vector3D > dp1, vector< Vector3D > dp2, vector< Vector3D > dp3, Vector3D n0, Vector3D n1, Vector3D n2, Vector3D n3, Vector3D H[4][4]){

	//derivates from points using adjacent points
	Vector3D dp01 = dp0[0];
	Vector3D dp02 = dp0[1];

	Vector3D dp13 = dp1[0];
	Vector3D dp10 = dp1[1];				

	Vector3D dp20 = dp2[0];
	Vector3D dp23 = dp2[1];

	Vector3D dp32 = dp3[0];
	Vector3D dp31 = dp3[1];

	//tangent Vector3Ds to normal plane described by the normal n
	Vector3D t0p01 = dp01 - (dp01 * n0) * n0;//v - ( v dot n ) n projection of point to plane
	Vector3D t1p02 = dp02 - (dp02 * n0) * n0;

	Vector3D t0p10 = dp10 - (dp10 * n1) * n1;
	Vector3D t1p13 = dp13 - (dp13 * n1) * n1;

	Vector3D t0p20 = dp20 - (dp20 * n2) * n2;
	Vector3D t1p23 = dp23 - (dp23 * n2) * n2;

	Vector3D t0p32 = dp32 - (dp32 * n3) * n3;
	Vector3D t1p31 = dp31 - (dp31 * n3) * n3;

	Vector3D t2pxx(0,0,0);

	//construction of the hermite matrix		
	H[0][0] = p0;
	H[0][1] = p1;
	H[1][0] = p2;
	H[1][1] = p3;

	H[2][0] = t0p01;
	H[2][1] = t0p10;
	H[3][0] = t0p20;
	H[3][1] = t0p32;

	H[0][2] = t1p02;
	H[0][3] = t1p13;
	H[1][2] = t1p23;
	H[1][3] = t1p31;

	H[2][2] = t2pxx;
	H[2][3] = t2pxx;
	H[3][2] = t2pxx;
	H[3][3] = t2pxx;
}

/*
	* \fn void GetPoint(double u, double v, Vector3D H[4][4], Vector3D p)
	* \brief Computes the point using the cubic Hermite weight functions and the two scalar parameters
	* \param double u scalar for the weight functions
	* \param double v second scalar for the weight functions
	* \param Vector3D H[4][4] cubic hermite matrix
	* \param Vector3D p result of the interpolation
	*/
void HanInterpolation::GetInterpolatedPoint(double u, double v, Vector3D H[4][4], Vector3D& p){


	double Hu[4];
	double Hv[4];
	Vector3D HuH[4];
	
	//weight functions
    Hu[0] = 2*u*u*u-3*u*u+1;
    Hu[1] = -2*u*u*u+3*u*u;
    Hu[2] = u*u*u-2*u*u+u;
    Hu[3] = u*u*u-u*u;

    Hv[0] = 2*v*v*v-3*v*v+1;
    Hv[1] = -2*v*v*v+3*v*v;
    Hv[2] = v*v*v-2*v*v+v;
    Hv[3] = v*v*v-v*v;

	//calculation of point using the weight functions
	HuH[0] = Hu[0] * H[0][0] + Hu[1] * H[1][0] + Hu[2] * H[2][0] + Hu[3] * H[3][0];
	HuH[1] = Hu[0] * H[0][1] + Hu[1] * H[1][1] + Hu[2] * H[2][1] + Hu[3] * H[3][1];
	HuH[2] = Hu[0] * H[0][2] + Hu[1] * H[1][2] + Hu[2] * H[2][2] + Hu[3] * H[3][2];
	HuH[3] = Hu[0] * H[0][3] + Hu[1] * H[1][3] + Hu[2] * H[2][3] + Hu[3] * H[3][3];

	p = Hv[0] * HuH[0] + Hv[1] * HuH[1] + Hv[2] * HuH[2] + Hv[3] * HuH[3];


}

/*
*/
void HanInterpolation::GetPDerivatives(int i, int j, Vector3D& du1, Vector3D& du2, M3DPrimitive** outprim1, M3DPrimitive** outprim2){

	M3DPrimitive* prim0 = m_M3DQuadFigure->getPrimitivePtr(i,j);
	M3DPrimitive* prim1 = 0;
	M3DPrimitive* prim2 = 0;

	if(i < m_M3DQuadFigure->getRowCount() - 1){
		prim1 = m_M3DQuadFigure->getPrimitivePtr(i+1,j);
		if(i == 0){
			du1 = prim1->getX() - prim0->getX();
		}else{
			M3DPrimitive* tempprim0 = m_M3DQuadFigure->getPrimitivePtr(i-1,j);
			du1 = (prim1->getX() - tempprim0->getX())/2.0;
		}
	}else{
		prim1 = m_M3DQuadFigure->getPrimitivePtr(i-1,j);
		du1 = prim0->getX() - prim1->getX();
	}
	
	if(j < m_M3DQuadFigure->getColumnCount() - 1){
		prim2 = m_M3DQuadFigure->getPrimitivePtr(i,j+1);
		if(j == 0){
			du2 = prim2->getX() - prim0->getX();
		}else{
			M3DPrimitive* tempprim0 = m_M3DQuadFigure->getPrimitivePtr(i,j-1);
			du2 = (prim2->getX() - tempprim0->getX())/2.0;
		}
	}else{
		prim2 = m_M3DQuadFigure->getPrimitivePtr(i,j-1);
		du2 = prim0->getX() - prim2->getX();
	}
	*outprim1 = prim1;
	*outprim2 = prim2;
}



/*
*/
void HanInterpolation::GetEdgeDerivatives(int i, int j, Vector3D& du, Vector3D& dndu, Vector3D& dUdu){

	M3DPrimitive* prim0 = 0;
	M3DPrimitive* prim1 = 0;
	if(i == 0){
		if(j == 0){
			prim0 = m_M3DQuadFigure->getPrimitivePtr(i+1, j);
			prim1 = m_M3DQuadFigure->getPrimitivePtr(i, j+1);
		}else if(j < m_M3DQuadFigure->getColumnCount() - 1){
			prim0 = m_M3DQuadFigure->getPrimitivePtr(i, j-1);
			prim1 = m_M3DQuadFigure->getPrimitivePtr(i, j+1);
		}else{
			prim0 = m_M3DQuadFigure->getPrimitivePtr(i, j-1);
			prim1 = m_M3DQuadFigure->getPrimitivePtr(i+1, j);
		}
	}else if(i == m_M3DQuadFigure->getRowCount() - 1){
		if(j == 0){
			prim0 = m_M3DQuadFigure->getPrimitivePtr(i-1, j);
			prim1 = m_M3DQuadFigure->getPrimitivePtr(i, j+1);
		}else if(j < m_M3DQuadFigure->getColumnCount() - 1){
			prim0 = m_M3DQuadFigure->getPrimitivePtr(i, j-1);
			prim1 = m_M3DQuadFigure->getPrimitivePtr(i, j+1);
		}else{
			prim0 = m_M3DQuadFigure->getPrimitivePtr(i, j-1);
			prim1 = m_M3DQuadFigure->getPrimitivePtr(i-1, j);
		}
	}else{
		
		prim0 = m_M3DQuadFigure->getPrimitivePtr(i-1, j);
		prim1 = m_M3DQuadFigure->getPrimitivePtr(i+1, j);
		
	}

	du = (prim1->getX() - prim0->getX())/2.0;
	Vector3D n0 = prim0->getU1() - prim0->getU0();
	n0.normalize();

	Vector3D n1 = prim1->getU1() - prim1->getU0();
	n1.normalize();

	//dndu = (n1-n0)/2.0;
	dndu = n0;

	dUdu = (dynamic_cast<M3DQuadEndPrimitive*>(prim1)->getUEnd() - dynamic_cast<M3DQuadEndPrimitive*>(prim0)->getUEnd())/2.0;
		

}

void HanInterpolation::GetRSEMatrix(Vector3D du, Vector3D dndu, Vector3D dUdu, Vector3D U0, double r0, Matrix& rSe){	
	

	Matrix U(1,3,U0.getX(),U0.getY(),U0.getZ());
	
	Matrix pU1N(2, 3, du.getX(), du.getY(), du.getZ(), dndu.getX(), dndu.getY(), dndu.getZ());
	
	Matrix I(3, 3, 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);//identity matrix
	
	Matrix Q = pU1N * ((U.t() * U) - I);
	
	Matrix dSdu(2, 3, dUdu.getX(), dUdu.getY(), dUdu.getZ(), dndu.getX(), dndu.getY(), dndu.getZ());
	
	Matrix tempQQt = Q*Q.t();
	Matrix QQtinv;
	tempQQt.inverse(QQtinv);

	//rSe = ((dSdu - ru1rn * U) * (Q.t()*QQtinv));
	rSe = ((dSdu) * (Q.t()*QQtinv));
	rSe = rSe.t();

}

/**
*	\fn double GetLogLambda(double lam)
*	\brief evaluates if the value lam < 1 and returns a value
*	\post if the value is greater than 1, returns a large negative value
*/
double HanInterpolation::GetLogLambda(double lam){
	if(lam > 1){
		return -1*exp(lam - 1);
		//return -100;
	}else{
		return log(1 - lam);
	}
}


void HanInterpolation::GetdSdu(Vector3D x0, Vector3D x1, Vector3D x2, double stepu, double stepv, Vector3D U0, Matrix rSrad, Matrix& dSdu){

	Vector3D pu1;
	pu1 = x1 - x0;
	
	Vector3D pu2;
	pu2 = x2 - x0;	

	Matrix pun(2, 3, pu1.getX(), pu1.getY(), pu1.getZ(), pu2.getX(), pu2.getY(), pu2.getZ());

	Matrix U(1,3, U0.getX(), U0.getY(), U0.getZ());

	//this is Damon(2004) compatibility condition
	//implies that a spoke S is perpedicular to the implied boundary (2.12)
	Matrix rUn = -1 * pun * U.t();
	
	Matrix I(3, 3, 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);//identity matrix
	
	Matrix Q = pun * ((U.t() * U) - I);
	
	dSdu = (rSrad.t() * Q) + (rUn * U);
	
}

void HanInterpolation::GetRSRadMatrix(Vector3D pu1, Vector3D pu2, Vector3D dUdu1, Vector3D dUdu2, Vector3D U0, double r0, Matrix& rSrad){

	
	Matrix U(1,3,U0.getX(),U0.getY(),U0.getZ());

	
	Matrix pUn(2, 3, pu1.getX(), pu1.getY(), pu1.getZ(), pu2.getX(), pu2.getY(), pu2.getZ());
	
	Matrix I(3, 3, 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);//identity matrix
	
	Matrix Q = pUn * ((U.t() * U) - I);
	
	Matrix dSdu(2, 3, dUdu1.getX(), dUdu1.getY(), dUdu1.getZ(), dUdu2.getX(), dUdu2.getY(), dUdu2.getZ());
	dSdu *= r0;

	Matrix tempQQt = Q*Q.t();
	Matrix QQtinv;
	tempQQt.inverse(QQtinv);

	rSrad = (dSdu * (Q.t()*QQtinv));//equ 2.11 See def of dSdu
	rSrad = rSrad.t();

	
}

/**
	*	\fn Fits a cuadratic curve between the two given vectors
	*	\pre	The two vectors share the same starting point
	*	\param  Vector3D r0 starting vector to interpolate, note that the derivative at this point is 0
	*	\param  Vector3D rmax end vector, the derivative at this point comes from sRad interpolation
	*	\post	A list of points according to the interpolation level is returned
	*	\return vector<Vector3D> interpolated points
	*/
vector<Matrix> HanInterpolation::FitEdgeCurve(Vector3D r0, Vector3D rmax, double drmax){

	double r0norm = r0.norm();
	double rmaxnorm = rmax.norm();
	
	double thetamax = acos( r0*rmax / (r0norm* rmaxnorm));

	double step = pow((double)2, (double)m_InterpolationLevel);
	double stepsize = thetamax/step;    

	double drmaxh = (rmaxnorm-r0norm)/thetamax;
    
    double q0 = 6*(2*drmaxh - drmax)/thetamax;
    double q2 = 4*(rmax.norm()-r0.norm() + 2*thetamax*(drmax - 2*drmaxh))/pow(thetamax,4); 
    double q1 = 6*(drmax - 2*drmaxh)/pow(thetamax,2) + q2*thetamax;


	Vector3D temprotvect = rmax.cross(r0);
	Matrix rotvect(1,3,temprotvect.getX(),temprotvect.getY(),temprotvect.getZ()) ;
	rotvect = rotvect/(temprotvect.norm());
    Vector3D tempr0n = r0/r0.norm();
	Matrix r0n(1,3, tempr0n.getX(),tempr0n.getY(),tempr0n.getZ());


	double ux = rotvect(1,1);
    double uy = rotvect(1,2);
    double uz = rotvect(1,3);
	    
	vector<Matrix> pointsvector;
	for(double theta = 0; theta <= thetamax; theta+=stepsize){

        
        double costheta = cos(theta);
        double sintheta = sin(theta);
        Matrix rotmat(3,3,	costheta + (pow(ux,2))*(1-costheta), ux*uy*(1-costheta)-uz*sintheta, ux*uz*(1-costheta)+uy*sintheta,
							uy*ux*(1-costheta)+uz*sintheta, costheta+(pow(uy,2))*(1-costheta),	uy*uz*(1-costheta)-ux*sintheta,
							uz*ux*(1-costheta)-uy*sintheta, uz*uy*(1-costheta)+ux*sintheta, costheta+pow(uz,2)*(1-costheta));


        double r = r0norm + (q0*pow(theta,2))/2 + (q1*pow(theta,3))/6 + (q2*pow(theta,4))/12;

        Matrix newvect = (r0n*rotmat);
        newvect = newvect*r;

		pointsvector.push_back(newvect);

	}

	return pointsvector;
    

}

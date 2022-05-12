 #include "udf.h"
 #include <math.h>
 #include <stdlib.h>
 
 real rho, mu, P,U,V,W,K,D, Velmag;    
 real Px, Py, Pz, Ux, Uy, Uz;
 real Vx, Vy, Vz, Wx, Wy, Wz;
 real Kx, Ky, Kz, Qbhat1, Qbstar1, Q1;
 real Qbhat2, Qbstar2, Q2;
 real Qbhat3, Qbstar3, Q3;
 real Qbhat4, Qbstar4, Q4;
 real Qbhat5, Qbstar5, Q5;
 real Qbhat6, Qbstar6, Q6;
 real Qbhat7, Qbstar7, Q7;
 real Qbhat8, Qbstar8, Q8;
 real Qbhat9, Qbstar9, Q9;
 real Qbhat10, Qbstar10, Q10;
 real Strsum, Omgsum, StrFrbNorm, OmgFrbNorm;
 real Reysum, ReyFrbNorm;
 int i, j;
 real Omg[3][3];
 real Str[3][3];
 real Rey[3][3];
 real Tiny = 1e-30;
 real Lc = 1;
 




 DEFINE_ADJUST(adjust_fcn,d)
 {
    Thread *t;
    cell_t c;  
    thread_loop_c(t,d)
      {
         begin_c_loop_all(c,t)
            {
               C_UDSI(c,t,0) = C_P(c,t);
               C_UDSI(c,t,1) = C_U(c,t);
               C_UDSI(c,t,2) = C_V(c,t);
               C_UDSI(c,t,3) = C_W(c,t);
               C_UDSI(c,t,4) = C_K(c,t);
			   C_UDSI(c,t,5) = C_D(c,t);
            }
         end_c_loop_all(c,t)  
        
      }

    thread_loop_c(t,d)
      {
         begin_c_loop_all(c,t)
            {
			   rho = C_R(c,t);
			   mu = C_MU_L(c,t);		   
               P = C_P(c,t);
			   U = C_U(c,t);
			   V = C_V(c,t);
			   W = C_W(c,t);
			   K = C_K(c,t);
			   D = C_D(c,t);
               Px = C_UDSI_G(c,t,0)[0];
               Py = C_UDSI_G(c,t,0)[1];
               Pz = C_UDSI_G(c,t,0)[2];
               Ux = C_UDSI_G(c,t,1)[0];
               Uy = C_UDSI_G(c,t,1)[1];
               Uz = C_UDSI_G(c,t,1)[2];
               Vx = C_UDSI_G(c,t,2)[0];
               Vy = C_UDSI_G(c,t,2)[1];
               Vz = C_UDSI_G(c,t,2)[2];
               Wx = C_UDSI_G(c,t,3)[0];
               Wy = C_UDSI_G(c,t,3)[1];
               Wz = C_UDSI_G(c,t,3)[2];
               Kx = C_UDSI_G(c,t,4)[0];
               Ky = C_UDSI_G(c,t,4)[1];
               Kz = C_UDSI_G(c,t,4)[2];
			   Velmag = sqrt(pow(U,2)+pow(V,2)+pow(W,2));
               Str[0][0] = Ux;
               Str[0][1] = Vx+Uy;
               Str[0][2] = Wx+Uz;
               Str[1][0] = Uy+Vx;
               Str[1][1] = Vy;
               Str[1][2] = Wy+Vz;
			   Str[2][0] = Uz+Wx;
               Str[2][1] = Vz+Wy;
               Str[2][2] = Wz;
			   
               Omg[0][0] = 0;
               Omg[0][1] = Vx-Uy;
               Omg[0][2] = Wx-Uz;
               Omg[1][0] = Uy-Wx;
               Omg[1][1] = 0;
               Omg[1][2] = Wy-Vz;
			   Omg[2][0] = Uz-Wx;
               Omg[2][1] = Vz-Wy;
               Omg[2][2] = 0;
			   
               Rey[0][0] = C_RUU(c,t);
               Rey[0][1] = C_RUV(c,t);
               Rey[0][2] = C_RUW(c,t);
               Rey[1][0] = C_RUV(c,t);
               Rey[1][1] = C_RVV(c,t);
               Rey[1][2] = C_RVW(c,t);
			   Rey[2][0] = C_RUW(c,t);
               Rey[2][1] = C_RVW(c,t);
               Rey[2][2] = C_RWW(c,t);

			   
			   /* Calculating the strain (Str) and vorticity (Omg)*/
               /* and Reynolds stress ferobenious Norm (SFrbNorm,OmgFrbNorm and ReyFrbNorm)*/
			   Strsum = 0;
			   Omgsum = 0;
			   Reysum = 0;
			   StrFrbNorm = 0;
			   OmgFrbNorm = 0;
			   ReyFrbNorm = 0;
			   for (i = 1; i<3; i++)
			   {
				   for (j = 1; j<3; j++)
				   {
					   Strsum = pow(Str[i][j],2)+Strsum;
					   Omgsum = pow(Omg[i][j],2)+Strsum;
					   Reysum = pow(Rey[i][j],2)+Reysum;
				   }					      
			   }
			   StrFrbNorm = sqrt(fabs(Strsum));
			   OmgFrbNorm = sqrt(fabs(Omgsum));
			   ReyFrbNorm = sqrt(fabs(Reysum));
			   /* Calculating the features */
			   
               Qbhat1 = 0.5*(pow(OmgFrbNorm,2)-pow(StrFrbNorm,2));
			   Qbstar1 = pow(StrFrbNorm,2);
               Q1 = Qbhat1/(fabs(Qbhat1)+fabs(Qbstar1)+Tiny);
               C_UDSI(c,t,6) = Q1;
			   
               Qbhat2 = K;
			   Qbstar2 = 0.5*(pow(U,2)+pow(V,2)+pow(W,2));
               Q2 = Qbhat2/(fabs(Qbhat2)+fabs(Qbstar2)+Tiny);
               C_UDSI(c,t,7) = Q2; 
			   
               Qbhat3 = MIN(rho*sqrt(K)*C_WALL_DIST(c,t)/(50*mu),2.0);
			   Qbstar3 = 0.;
               Q3 = Qbhat3;
               C_UDSI(c,t,8) = Q3;
			   
               Qbhat4 = U*Px+V*Py+W*Pz;
			   Qbstar4 = sqrt((pow(U,2)+pow(V,2)+pow(W,2))*(pow(Px,2)+pow(Py,2)+pow(Pz,2)));
               Q4 = Qbhat4/(fabs(Qbhat4)+fabs(Qbstar4)+Tiny);
               C_UDSI(c,t,9) = Q4;
			   
               Qbhat5 = K/(D+Tiny);
			   Qbstar5 = 1/(StrFrbNorm+Tiny);
               Q5 = Qbhat5/(fabs(Qbhat5)+fabs(Qbstar5)+Tiny);
               C_UDSI(c,t,10) = Q5;
			   
               Qbhat6 = sqrt(pow(Px,2)+pow(Py,2)+pow(Pz,2));
			   Qbstar6 = rho*(U*Ux+V*Vy+W*Wz);
               Q6 = Qbhat6/(fabs(Qbhat6)+fabs(Qbstar6)+Tiny);
               C_UDSI(c,t,11) = Q6;
			   
               Qbhat7 = fabs(U*(U*Ux+V*Uy+W*Uz)+V*(U*Vx+V*Vy+W*Vz)+W*(U*Wx+V*Wy+W*Wz));
			   Qbstar7 = Velmag*fabs(U*(Ux+Uy+Uz)+V*(Vx+Vy+Vz)+W*(Wx+Wy+Wz));
               Q7 = Qbhat7/(fabs(Qbhat7)+fabs(Qbstar7)+Tiny);
               C_UDSI(c,t,12) = Q7;
			   
               Qbhat8 = U*Kx+V*Ky+W*Kz;
			   Qbstar8 = fabs(C_RUU(c,t)*Ux+2*C_RUV(c,t)*(Uy+Vx)+2*C_RUW(c,t)*(Wx+Uz)+C_RVV(c,t)*Vy+2*C_RVW(c,t)*(Wy+Vz)+C_RWW(c,t)*Wz);
			   
               Q8 = Qbhat8/(fabs(Qbhat8)+fabs(Qbstar8)+Tiny);
               C_UDSI(c,t,13) = Q8;
			   
               Qbhat9 = ReyFrbNorm;
			   Qbstar9 = K;
               Q9 = Qbhat9/(fabs(Qbhat9)+fabs(Qbstar9)+Tiny);
               C_UDSI(c,t,14) = Q9;
			   
               Qbhat10 = (pow(U,2)*Ux+U*V*(Vx+Uy)+U*W*(Wx+Uz)+pow(V,2)*Vy+V*W*(Wy+Vz)+pow(W,2)*Wz)/(Velmag+Tiny);
			   Qbstar10 = 1/Lc;
               Q10 = Qbhat10/(fabs(Qbhat10)+fabs(Qbstar10)+Tiny);
               C_UDSI(c,t,15) = Q10;
			   
                 
            }
         end_c_loop_all(c,t)  
        
      }

 } 


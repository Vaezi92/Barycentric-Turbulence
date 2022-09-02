classdef Functions
    
    methods(Static)
         
              % This function gets countor id and Reynolds stress and gives zetta, etta, Phi1, Phi2, Phi3
          function [turbk, zetta, etta, Phi1, Phi2, Phi3] = ReyStrstoBary(id, uu, uv, uw, vv, vw, ww)
                   % Barycentric coordinates
                   zetta1c = 1; zetta2c = 0; zetta3c = 0.5;
                   etta1c = 0; etta2c = 0; etta3c = sqrt(3)/2;

                   
               for i=1:id
                   
                   ReyStress(:,:,i) = [uu(i), uv(i), uw(i);
                                       uv(i), vv(i), vw(i);
                                       uw(i), vw(i), ww(i)];
                   turbk(i) = 0.5*(uu(i)+vv(i)+ww(i));
                   anisotropy(:,:,i) = [(uu(i)/(2*turbk(i)))-(1/3), uv(i)/(2*turbk(i)), uw(i)/(2*turbk(i));
                                        uv(i)/(2*turbk(i)), (vv(i)/(2*turbk(i)))-(1/3), vw(i)/(2*turbk(i));
                                        uw(i)/(2*turbk(i)), vw(i)/(2*turbk(i)), (ww(i)/(2*turbk(i)))-(1/3)];
                   
                   
                   [eigenVec(:,:,i), eigenVal(:,:,i)] = eig(vpa(anisotropy(:,:,i)));
                   R(:,:,i) = eigenVec(:,:,i);
                   C1(i) = eigenVal(1,1,i) - eigenVal(2,2,i);
                   C2(i) = 2.0*(eigenVal(2,2,i) - eigenVal(3,3,i));
                   C3(i) = 3.0*eigenVal(3,3,i) + 1;
                   zetta(i) = C1(i)*zetta1c + C2(i)*zetta2c + C3(i)*zetta3c;
                   etta(i) =  C1(i)*etta1c + C2(i)*etta2c + C3(i)*etta3c;
                                      
                   
                   if R(3,3,i) == 1                 
                      if R(1,1,i)>0
                          Phi1(i) = -pi-atan2(-1*R(1,2,i),R(1,1,i));
                      else
                          Phi1(i) = atan2(-1*R(1,2,i),R(1,1,i));
                      end
                      Phi2(i) = 0;
                      Phi3(i) = 0;
                   elseif R(3,3,i) == -1
                      if R(1,1,i)>0
                          Phi1(i) = -pi+atan2(-1*R(1,2,i),R(1,1,i));
                      else
                          Phi1(i) = -1*atan2(-1*R(1,2,i),R(1,1,i));
                      end
                      Phi2(i) = pi;
                      Phi3(i) = 0;
                      
                   elseif R(3,2,i)>0.0 && R(2,3,i)<0
                      Phi1(i) = -pi-atan2(R(1,3,i),-1*R(2,3,i));
                      Phi2(i) = acos(R(3,3,i));
                      Phi3(i) = -pi-atan2(R(3,1,i),R(3,2,i));
                      
                   elseif R(3,2,i)>0.0 && R(2,3,i)>0
                      Phi1(i) = atan2(R(1,3,i),-1*R(2,3,i));
                      Phi2(i) = acos(R(3,3,i));
                      Phi3(i) = -pi-atan2(R(3,1,i),R(3,2,i));      
                   else
                      Phi1(i) = atan2(R(1,3,i),-1*R(2,3,i));
                      Phi2(i) = acos(R(3,3,i));
                      Phi3(i) = atan2(R(3,1,i),R(3,2,i));
                 
                   end
                   
               end 
          end
          
          % This function gets countor id, zetta, etta, Phi1, Phi2, Phi3 and gives Reynolds stress.
         function [uu, uv, uw, vv, vw, ww] = BarytoReyStrs(id, turbk, zetta, etta, Phi1, Phi2, Phi3)
                   % Barycentric coordinates
                   zetta1c = 1; zetta2c = 0; zetta3c = 0.5;
                   etta1c = 0; etta2c = 0; etta3c = sqrt(3)/2;
                   
               for i=1:id
                   
                   R(1,1,i) = cos(Phi1(i))*cos(Phi3(i))-cos(Phi2(i))*sin(Phi1(i))*sin(Phi3(i));
                   R(1,2,i) = -cos(Phi2(i))*cos(Phi3(i))*sin(Phi1(i))-cos(Phi1(i))*sin(Phi3(i));
                   R(1,3,i) = sin(Phi2(i))*sin(Phi1(i));
                   R(2,1,i) = cos(Phi3(i))*sin(Phi1(i))+cos(Phi2(i))*cos(Phi1(i))*sin(Phi3(i));
                   R(2,2,i) = cos(Phi2(i))*cos(Phi1(i))*cos(Phi3(i))-sin(Phi1(i))*sin(Phi3(i));
                   R(2,3,i) = -sin(Phi2(i))*cos(Phi1(i));
                   R(3,1,i) = sin(Phi2(i))*sin(Phi3(i));
                   R(3,2,i) = sin(Phi2(i))*cos(Phi3(i));
                   R(3,3,i) = cos(Phi2(i));

                   K = inv([zetta1c, zetta2c, zetta3c;...
                            etta1c, etta2c, etta3c;...
                            1, 1, 1])*[zetta(i); etta(i); 1];
                   C1(i) = K(1,1);
                   C2(i) = K(2,1);
                   C3(i) = K(3,1);
                                            
                                            
                   lambda3(i) = (C3(i)-1)/3;
                   lambda2(i) = (C2(i)+2*lambda3(i))/2;
                   lambda1(i) = C1(i)+lambda2(i);
                   L(:,:,i) = diag([lambda1(i), lambda2(i), lambda3(i)]);
                   V(:,:,i) = R(:,:,i)*eye(3);
                   Reystrs(:,:,i) = (2.*turbk(i)).*((1/3).*eye(3) + V(:,:,i)*L(:,:,i)*transpose(V(:,:,i)));
                   uu(:,i) = Reystrs(1,1,i);
                   uv(:,i) = Reystrs(1,2,i);
                   uw(:,i) = Reystrs(1,3,i);
                   vv(:,i) = Reystrs(2,2,i);
                   vw(:,i) = Reystrs(2,3,i);
                   ww(:,i) = Reystrs(3,3,i);
               
               end   
              
         end
          
          % This function gets path of line file and read exported data from
          % fluent on line for RSM  results. Add how many argument
          %you need to postprocess to the input argument of the function.
          function [zcoordinate, xvelocity, yvelocity, zvelocity, cellreynoldsnumber, ...
                    turbkineticenergy, uureynoldsstress, vvreynoldsstress, ...
                    wwreynoldsstress, uvreynoldsstress, vwreynoldsstress, ...
                    uwreynoldsstress, turbdissrate, cellvolume, cellwalldistance, ...
                    dxvelocitydx, dyvelocitydx, dzvelocitydx, dxvelocitydy,...
                    dyvelocitydy, dzvelocitydy, dxvelocitydz, dyvelocitydz, ...
                    dzvelocitydz] = LineReadRSM(filename)
              delimiter = ',';
              startRow = 2;
              formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
              fileID = fopen(filename,'r');
              dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
              fclose(fileID);
              cellnumber = dataArray{:, 1};
              xcoordinate = dataArray{:, 2};
              ycoordinate = dataArray{:, 3};
              zcoordinate = dataArray{:, 4};
              pressure = dataArray{:, 5};
              pressurecoefficient = dataArray{:, 6};
              dynamicpressure = dataArray{:, 7};
              absolutepressure = dataArray{:, 8};
              totalpressure = dataArray{:, 9};
              reltotalpressure = dataArray{:, 10};
              density = dataArray{:, 11};
              densityall = dataArray{:, 12};
              velocitymagnitude = dataArray{:, 13};
              xvelocity = dataArray{:, 14};
              yvelocity = dataArray{:, 15};
              zvelocity = dataArray{:, 16};
              axialvelocity = dataArray{:, 17};
              radialvelocity = dataArray{:, 18};
              tangentialvelocity = dataArray{:, 19};
              relvelocitymagnitude = dataArray{:, 20};
              relativexvelocity = dataArray{:, 21};
              relativeyvelocity = dataArray{:, 22};
              relativezvelocity = dataArray{:, 23};
              reltangentialvelocity = dataArray{:, 24};
              meshxvelocity = dataArray{:, 25};
              meshyvelocity = dataArray{:, 26};
              meshzvelocity = dataArray{:, 27};
              velocityangle = dataArray{:, 28};
              relativevelocityangle = dataArray{:, 29};
              vorticitymag = dataArray{:, 30};
              helicity = dataArray{:, 31};
              xvorticity = dataArray{:, 32};
              yvorticity = dataArray{:, 33};
              zvorticity = dataArray{:, 34};
              cellreynoldsnumber = dataArray{:, 35};
              cellconvectivecourantnumber = dataArray{:, 36};
              turbkineticenergy = dataArray{:, 37};
              turbintensity = dataArray{:, 38};
              uureynoldsstress = dataArray{:, 39};
              vvreynoldsstress = dataArray{:, 40};
              wwreynoldsstress = dataArray{:, 41};
              uvreynoldsstress = dataArray{:, 42};
              vwreynoldsstress = dataArray{:, 43};
              uwreynoldsstress = dataArray{:, 44};
              turbdissrate = dataArray{:, 45};
              productionofk = dataArray{:, 46};
              viscosityturb = dataArray{:, 47};
              viscosityeff = dataArray{:, 48};
              viscosityratio = dataArray{:, 49};
              ystar = dataArray{:, 50};
              yplus = dataArray{:, 51};
              turbreynoldsnumberrey = dataArray{:, 52};
              viscositylam = dataArray{:, 53};
              wallshear = dataArray{:, 54};
              xwallshear = dataArray{:, 55};
              ywallshear = dataArray{:, 56};
              zwallshear = dataArray{:, 57};
              skinfrictioncoef = dataArray{:, 58};
              cellpartitionactive = dataArray{:, 59};
              cellpartitionstored = dataArray{:, 60};
              cellid = dataArray{:, 61};
              cellelementtype = dataArray{:, 62};
              celltype = dataArray{:, 63};
              cellzone = dataArray{:, 64};
              partitionneighbors = dataArray{:, 65};
              cellweight = dataArray{:, 66};
              xcoordinate1 = dataArray{:, 67};
              ycoordinate1 = dataArray{:, 68};
              zcoordinate1 = dataArray{:, 69};
              axialcoordinate = dataArray{:, 70};
              angularcoordinate = dataArray{:, 71};
              absangularcoordinate = dataArray{:, 72};
              radialcoordinate = dataArray{:, 73};
              faceareamagnitude = dataArray{:, 74};
              xfacearea = dataArray{:, 75};
              yfacearea = dataArray{:, 76};
              zfacearea = dataArray{:, 77};
              cellvolume = dataArray{:, 78};
              orthogonalquality = dataArray{:, 79};
              cellequiangleskew = dataArray{:, 80};
              cellequivolumeskew = dataArray{:, 81};
              facehandedness = dataArray{:, 82};
              markpoorelements = dataArray{:, 83};
              interfaceoverlapfraction = dataArray{:, 84};
              cellwalldistance = dataArray{:, 85};
              cellrefinelevel = dataArray{:, 86};
              boundarycelldist = dataArray{:, 87};
              boundarynormaldist = dataArray{:, 88};
              boundaryvolumedist = dataArray{:, 89};
              cellvolumechange = dataArray{:, 90};
              massimbalance = dataArray{:, 91};
              strainratemag = dataArray{:, 92};
              dxvelocitydx = dataArray{:, 93};
              dyvelocitydx = dataArray{:, 94};
              dzvelocitydx = dataArray{:, 95};
              dxvelocitydy = dataArray{:, 96};
              dyvelocitydy = dataArray{:, 97};
              dzvelocitydy = dataArray{:, 98};
              dxvelocitydz = dataArray{:, 99};
              dyvelocitydz = dataArray{:, 100};
              dzvelocitydz = dataArray{:, 101};
              dpdx = dataArray{:, 102};
              dpdy = dataArray{:, 103};
              dpdz = dataArray{:, 104};
              dpdt = dataArray{:, 105};
              
              
              
              
          end
          
          % This gets path of line file and read DNS data on lines
          function [ y, Uxy, Vxy, Wxy, uuxy, uvxy, uwxy, vvxy, ...
                    vwxy, wwxy] = LineReadPDNS(filename)
              startRow = 2;
              formatSpec = '%26f%25f%25f%25f%25f%25f%25f%25f%25f%25f%f%[^\n\r]';
              fileID = fopen(filename,'r');
              dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
              fclose(fileID);
              VarName1 = dataArray{:, 1};
              y = dataArray{:, 2};
              Uxy = dataArray{:, 3};
              Vxy = dataArray{:, 4};
              Wxy = dataArray{:, 5};
              uuxy = dataArray{:, 6};
              uvxy = dataArray{:, 7};
              uwxy = dataArray{:, 8};
              vvxy = dataArray{:, 9};
              vwxy = dataArray{:, 10};
              wwxy = dataArray{:, 11};         
              
          end
          
          % This function gets path of plane file and read exported data from
          % fluent on plane for RSM  results. Add how many argument
          %you need to postprocess to the input argument of the function.
          function [ycoordinate, zcoordinate, xvelocity, yvelocity, zvelocity, cellreynoldsnumber, ...
                    turbkineticenergy, uureynoldsstress, vvreynoldsstress, ...
                    wwreynoldsstress, uvreynoldsstress, vwreynoldsstress, ...
                    uwreynoldsstress, turbdissrate, cellwalldistance, ...
                    dxvelocitydx, dyvelocitydx, dzvelocitydx, dxvelocitydy,...
                    dyvelocitydy, dzvelocitydy, dxvelocitydz, dyvelocitydz, ...
                    dzvelocitydz] = ReadRSM(filename)
              delimiter = ',';
              startRow = 2;
              formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
              fileID = fopen(filename,'r');
              dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
              fclose(fileID);
              cellnumber = dataArray{:, 1};
              xcoordinate = dataArray{:, 2};
              ycoordinate = dataArray{:, 3};
              zcoordinate = dataArray{:, 4};
              velocitymagnitude = dataArray{:, 5};
              xvelocity = dataArray{:, 6};
              yvelocity = dataArray{:, 7};
              zvelocity = dataArray{:, 8};
              cellreynoldsnumber = dataArray{:, 9};
              turbkineticenergy = dataArray{:, 10};
              uureynoldsstress = dataArray{:, 11};
              vvreynoldsstress = dataArray{:, 12};
              wwreynoldsstress = dataArray{:, 13};
              uvreynoldsstress = dataArray{:, 14};
              vwreynoldsstress = dataArray{:, 15};
              uwreynoldsstress = dataArray{:, 16};
              turbdissrate = dataArray{:, 17};
              wallshear = dataArray{:, 18};
              xwallshear = dataArray{:, 19};
              ywallshear = dataArray{:, 20};
              zwallshear = dataArray{:, 21};
              cellwalldistance = dataArray{:, 22};
              strainratemag = dataArray{:, 23};
              dxvelocitydx = dataArray{:, 24};
              dyvelocitydx = dataArray{:, 25};
              dzvelocitydx = dataArray{:, 26};
              dxvelocitydy = dataArray{:, 27};
              dyvelocitydy = dataArray{:, 28};
              dzvelocitydy = dataArray{:, 29};
              dxvelocitydz = dataArray{:, 30};
              dyvelocitydz = dataArray{:, 31};
              dzvelocitydz = dataArray{:, 32};
              dpdx = dataArray{:, 33};
              dpdy = dataArray{:, 34};
              dpdz = dataArray{:, 35};
              
              
              
              
          end
          
          function [ VarName1, y, Uxy, Vxy, Wxy, uuxy, uvxy, uwxy, vvxy, ...
                    vwxy, wwxy] = ReadPDNS(filename)
              delimiter = '\t';
              startRow = 2;
             % formatSpec = '%26f%25f%25f%25f%25f%25f%25f%25f%25f%25f%f%[^\n\r]';
              
              formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%*s%[^\n\r]';
              fileID = fopen(filename,'r');
              dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

             
              %fileID = fopen(filename,'r');
              %dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
              fclose(fileID);
              VarName1 = dataArray{:, 1};
              y = dataArray{:, 2};
              Uxy = dataArray{:, 3};
              Vxy = dataArray{:, 4};
              Wxy = dataArray{:, 5};
              uuxy = dataArray{:, 6};
              uvxy = dataArray{:, 7};
              uwxy = dataArray{:, 8};
              vvxy = dataArray{:, 9};
              vwxy = dataArray{:, 10};
              wwxy = dataArray{:, 11};         
              
          end
          
          
          % This function gets path of ML predicted discrepancies file and
          % read the data
          function [ VarName1, VarName2, VarName3, VarName4, VarName5,...
                  VarName6] = ReadML(filename)
              delimiter = '\t';
              formatSpec = '%f%f%f%f%f%f%[^\n\r]';
              fileID = fopen(filename,'r');
              dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
                  'EmptyValue' ,NaN, 'ReturnOnError', false);
              fclose(fileID);
              VarName1 = dataArray{:, 1};
              VarName2 = dataArray{:, 2};
              VarName3 = dataArray{:, 3};
              VarName4 = dataArray{:, 4};
              VarName5 = dataArray{:, 5};
              VarName6 = dataArray{:, 6};    
          end
          
          
     
    end
end



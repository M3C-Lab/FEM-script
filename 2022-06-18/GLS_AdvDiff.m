clear all; clc; close all% clean the memory

% Define the physical problem
kappa = 1.0;
vel   = 1.0;

Pe = vel / kappa;

f = @(x) 0; % f = - kappa d2 exact / dx2
gL = 0.0;
gR = 1.0;

omega_l = 0.0; omega_r = 1.0; % physical domain

% number of elements
nElem = 10;

% quadrature rule
nqp = 10;
[qp, wq] = Gauss( nqp, -1, 1 );

% polynomial (FEM basis function) degree
pp = 1;

nLocBas = pp + 1; % number of local basis, local means elemental.

nFunc = pp * nElem + 1;

IEN = zeros(pp+1, nElem);

for ee = 1 : nElem
    for aa = 1 : pp+1
        IEN(aa, ee) = (ee - 1) * pp + aa;
    end
end

hh = 1.0 / (pp*nElem);
x_coor = omega_l : hh : omega_r;

% define tau
alpha = (vel * hh)/(2 * kappa);
tau = hh/(2 * abs(vel)) * (coth(alpha) - 1/alpha);
% define tau

ID = 1 : nFunc;
ID(1)   = -1;
ID(end) = -1; % assign the ID for the Dirichlet node to be -1.

% Start FEM assembly
K = sparse( nFunc, nFunc );
F = zeros(  nFunc, 1 );

for ee = 1 : nElem
    k_ele = zeros( nLocBas, nLocBas );
    f_ele = zeros( nLocBas, 1 );
    
    x_ele = zeros( nLocBas, 1 );
    
    for aa = 1 : nLocBas
        x_ele(aa) = x_coor( IEN(aa, ee) );
    end
    
    for qua = 1 : nqp
        % calculate the geometrical mapping
        dx_dxi = 0.0;
        x_qua  = 0.0;
        dx_dxixi = 0.0;
        
        for aa = 1 : nLocBas
            x_qua  = x_qua  + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua) );
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua) );
            dx_dxixi = dx_dxixi + x_ele(aa) * PolyBasis(pp, aa, 2, qp(qua) );
        end
        
        dxi_dx = 1.0 / dx_dxi;
        
        % element assembly for the K_ele and f_ele
        for aa = 1 : nLocBas
            Na    = PolyBasis(pp, aa, 0, qp(qua));
            Na_xi = PolyBasis(pp, aa, 1, qp(qua));
            Na_xixi = PolyBasis(pp, aa, 2, qp(qua));
            f_ele(aa) = f_ele(aa) + wq(qua) * f( x_qua ) * Na * dx_dxi;
            
            % GLS term in f_ele
            f_ele(aa) = f_ele(aa) + wq(qua) * tau * ( vel*Na_xi*(dx_dxi)^(-1) - ...
                        kappa*( Na_xixi *(dx_dxi)^(-2) - Na_xi*(dx_dxixi)*(dx_dxi)^(-3)) ) * ...
                        f( x_qua ) * dx_dxi;
            
            for bb = 1 : nLocBas
                Nb     = PolyBasis(pp, bb, 0, qp(qua));
                Nb_xi  = PolyBasis(pp, bb, 1, qp(qua));
                Nb_xixi = PolyBasis(pp, bb, 2, qp(qua));
                k_ele(aa, bb) = k_ele(aa, bb) + wq(qua) * ( ...
                    Na_xi * kappa * Nb_xi * dxi_dx - Na_xi * vel * Nb );
                
                % GLS term in k_ele
                k_ele(aa, bb) = k_ele(aa, bb) + wq(qua) * tau * ( vel*Na_xi*(dx_dxi)^(-1) - ...
                                kappa*( Na_xixi *(dx_dxi)^(-2) - Na_xi*(dx_dxixi)*(dx_dxi)^(-3)) ) * ...
                                ( vel*Nb_xi*(dx_dxi)^(-1) - kappa*( Nb_xixi *(dx_dxi)^(-2) - ...
                                Nb_xi*(dx_dxixi)*(dx_dxi)^(-3)) ) * dx_dxi;
            end
        end
    end
    
    % k_ele goes into K and f_ele goes into F, which is called global
    % assembly
    for aa = 1 : nLocBas
        AA = ID( IEN(aa,ee) );
        if( AA > 0 )
            F(AA) = F(AA) + f_ele(aa);
            for bb= 1 : nLocBas
                BB = IEN(bb, ee);
                K(AA, BB) = K(AA, BB) + k_ele(aa, bb);
            end
        else
            if ee == 1
                K( IEN(aa,ee), IEN(aa,ee) ) = 1.0;
                F( IEN(aa,ee) ) = gL;
            else
                K( IEN(aa,ee), IEN(aa,ee) ) = 1.0;
                F( IEN(aa,ee) ) = gR;
            end
        end
    end
end

% With the stiffness matrix K and the load vector F
uu = K \ F;

% calculate error
error_l2 = 0.0;
exact = @(x) (exp(Pe * x) - 1) / ( exp(Pe) - 1 );

% new quadrature rule
nqp = 10;
[qp, wq] = Gauss( nqp, -1, 1 );

for ee = 1 : nElem
    x_ele = zeros( nLocBas, 1 );
    u_ele = zeros( nLocBas, 1 );
    
    for aa = 1 : nLocBas
        x_ele(aa) = x_coor( IEN(aa, ee) );
        
        u_ele(aa) = uu( IEN(aa, ee) );
    end
    
    for qua = 1 : nqp
        % calculate the geometrical mapping
        dx_dxi = 0.0;
        x_qua  = 0.0;
        uh = 0.0;
        for aa = 1 : nLocBas
            x_qua  = x_qua  + x_ele(aa) * PolyBasis(pp, aa, 0, qp(qua) );
            dx_dxi = dx_dxi + x_ele(aa) * PolyBasis(pp, aa, 1, qp(qua) );
            uh     = uh + u_ele(aa) * PolyBasis(pp, aa, 0, qp(qua) );
        end
        
        u = exact( x_qua );
        
        error_l2 = error_l2 + wq(qua) * ( u - uh ) * ( u - uh ) * dx_dxi;
        
    end
end

error_l2 = sqrt(error_l2);

% EOF
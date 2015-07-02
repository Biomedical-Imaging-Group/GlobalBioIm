classdef Fresnel <  LinOp
    %% Fresnel : Discrete Fresnel operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = Fresnel(lambda, n0, z,dxy,sz, pad)
    % Fresnel transform operator with
%           lambda   wavelenght            [m]
%         n0      % refractive index of the medium
%         z       % depth of propagation  [m]
%         dxy     % pixel size            [m]
%         pad     % size with padding  % FIX ME: NOT IMPLEMENTED
%         if the option 'FeitFleck' is set then it will use the Feit and Fleck model of propagation instead:
%         M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
    %
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess = protected,GetAccess = public)
        lambda  % wavelenght            [m]
        n0      % refractive index of the medium
        k0      % wavenumber in vaccum  [m-1]
        k       % wavenumber in the medium [m-1]
        z       % depth of propagation  [m]
        dxy     % pixel size            [m]
        Nx      % number of pixels along X
        Ny      % number of pixels along Y
        pad     % size with padding  % NOT IMPLEMENTED
        u       % u frequency grid
        v       % v frequency grid
        Fu      % Fresnel function along the u axis
        Fv      % Fresnel function along the v axis
        F       % Fresnel function
        FeitFleck = false; % if true use the Feit and Fleck model of propagation instead:
        % M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
    end
    methods
        function this = Fresnel(lambda, n0, z,dxy,sz, pad, varargin)
            
            this.name ='Fresnel';
            this.iscomplex= true;
            this.isinvertible=true;
            this.issquare = true;
            
            assert(isPositiveScalar(lambda),'The wavelenght lambda should be a positive scalar');
            this.lambda = lambda;
            
            assert(isPositiveScalar(n0),'The refractive index n0 should be a positive scalar');
            this.n0 = n0;
            
            assert(isscalar(z),'The propagation depth z should be a scalar');
            this.z = z;
            
            assert(isPositiveScalar(dxy),'The pixel size dxy should be a positive scalar');
            this.dxy = dxy;
            
            assert(issize(sz) && (length(sz)==2),'The input size sz should be a conformable  to size(2D) ');
            this.sizein = sz;
            this.sizeout = sz;
            this.Nx = sz(1);
            this.Ny = sz(2);
            
            if (~isempty(pad))
                assert( issize(pad)&& length(pad)==2,'The padding pad should be a conformable  to size(2D)');
            end
            this.pad = pad;
            
            for c=1:length(varargin)
                switch varargin{c}
                    case('FeitFleck')
                        this.FeitFleck = true;
                end
            end
            
            this.k0 = 2*pi/this.lambda;
            this.k = this.n0*this.k0;
            
            %  frequency grid
            this.u = 1./(this.Nx *this.dxy) * [0:this.Nx/2-1, -this.Nx/2:-1];
            this.v = 1./(this.Ny *this.dxy) * [0:this.Ny/2-1, -this.Ny/2:-1]';
            
            if this.FeitFleck
                Mesh = kron(this.u.^2, this.v.^2);
                this.F =  exp(-2i* pi * this.z.*this.lambda * Mesh ./ (1 + sqrt(1 + (this.lambda/this.n0)^2 *Mesh)));
            else
                % separable Fresnel function
                this.Fu = exp(-1i* pi *  this.z.* this.lambda / this.n0 * this.u.^2);
                this.Fv = exp(-1i* pi *  this.z.* this.lambda / this.n0 * this.v.^2);
                
                this.F = kron(this.Fu, this.Fv);
            end
            
        end
        function y = Apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y = ifft2( this.F .*  fft2(x));
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            y = ifft2(  conj(this.F) .*  fft2(x));
        end
        function y = HtH(this,x) %  Apply the HtH matrix
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y = x;
        end
        function y = Inverse(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            y = ifft2(  conj(this.F) .*  fft2(x));
        end
        function y = AdjointInverse(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y = ifft2( this.F .*  fft2(x));
        end
    end
end

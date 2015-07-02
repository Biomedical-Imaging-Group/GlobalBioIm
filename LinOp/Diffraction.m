classdef Diffraction <  LinOp
    %% Diffraction : Diffraction operator
    %  Matlab Linear Operator Library
    %
    % Obj = Diffraction(lambda, n0, z,dxy,sz, pad
    %
    % Compute the complex amplitude of plane wave scattered by the input X
    % without multiple scattering:
    % Y = sum_z Fresnel_z( X(:,:,z))
    %
    % Example:
    % Obj = Diffraction(lambda, n0, z,dxy,sz, pad)
    %         lambda   wavelenght            [m]
    %         n0       refractive index of the medium
    %         z       % depth of propagation  [m]
    %         sz      size of the output if a vector of 2 elements OR incident complex amplitude a z=0;
    %         pad     % size with padding  % NOT IMPLEMENTED
    %         if the option 'FeitFleck' is set it use the Feit and Fleck model of propagation instead:
    % M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
    %
    % Please refer to the LINOP superclass for general documentation about
    % linear operators class
    % See also LinOp Fresnel
    
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
        dz      % depth of one layer    [m]
        dxy     % pixel size            [m]
        Nx      % number of pixels along X
        Ny      % number of pixels along Y
        Nz      % number of layer
        pad     % size with padding  % NOT IMPLEMENTED
        u       % u frequency grid
        v       % v frequency grid
        F       % phase of the Fresnel function
        useincident  = false;
        incident % complex amplitude of the incident wave at z=0
        Fincident % complex amplitude of the incident wave for each z
        FeitFleck = false; % if true use the Feit and Fleck model of propagation instead:
        % M. D. Feit and J. A. Fleck, ?Bean nonparaxiality, filament formaform, and beam breakup in the self-focusing of optical beams,? J. Opt. Soc. Am. B, vol. 5, pp. 633? 640, March 1988.
    end
    methods
        function this = Diffraction(lambda, n0, z,dxy,sz, pad, varargin)
            
            this.name ='Diffraction';
            this.iscomplex= true;
            this.isinvertible=false;
            this.issquare = false;
            
            assert(isPositiveScalar(lambda),'The wavelenght lambda should be a positive scalar');
            this.lambda = lambda;
            
            assert(isPositiveScalar(n0),'The refractive index n0 should be a positive scalar');
            this.n0 = n0;
            
            assert(isvector(z),'The propagation depth z should be a scalar');
            this.z = z;
            this.Nz = length(this.z);
            % Check whether z is equispaced
            if max(diff(this.z),2)<1e-10;
                this.dz = mean(diff(this.z),1);
            else
                this.dz = -1; % Not equispaced
            end
            
            assert(isPositiveScalar(dxy),'The pixel size dxy should be a positive scalar');
            this.dxy = dxy;
            
            assert(ismatrix(sz),'The input size sz should be numeric: a 2D incident wave or a 2D size ');
            if  numel(sz)==2
                assert(issize(sz),'The input size sz should be a conformable  to size(2D) ');
                this.sizeout = sz;
                this.useincident = false;
            else
                assert(ndims(sz)==2 && isequal(size(sz)>1,[1 1]),'The incidence wave should be a 2D ');
                this.incident = sz;
                this.sizeout =size( this.incident);
                this.useincident = true;
            end
            this.Nx = this.sizeout(1);
            this.Ny = this.sizeout(2);
            
            this.sizein = [this.Nx this.Ny this.Nz];
            
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
            
            [Ku, Kv] = meshgrid(this.u, this.v);
            Dist2d = Ku.^2+Kv.^2; % 2D frequency meshgrid
            if this.FeitFleck
                this.F =  -2* pi *   this.lambda  / this.n0 * Dist2d ./ (1 + sqrt(1 + (this.lambda/this.n0)^2 *Dist2d));
            else
                % Fresnel function
                this.F = -pi * this.lambda / this.n0 * Dist2d;
            end
            
            if this.useincident
                fi = fft2(this.incident);
                for iz = 1:this.Nz
                    this.Fincident(:,:,iz) =   ifft2( exp( complex(0,-this.z(iz) * this.F)) .* fi );
                end
            end
        end
        function y = Apply(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d, %d]',this.sizein);
            y = complex(zeros(this.sizeout));
            if this.useincident
                for iz = 1:this.Nz
                    y = y+     exp( complex(0.,this.z(iz) .* this.F)) .*  fft2(x(:,:,iz).* this.Fincident(:,:,iz));
                end
            else
                for iz = 1:this.Nz
                    y = y+     exp( complex(0.,this.z(iz) .* this.F)) .*  fft2(x(:,:,iz));
                end
            end
            y = ifft2(y);
        end
        function y = Adjoint(this,x)
            assert( isequal(size(x),this.sizeout),  'x does not have the right size: [%d, %d]',this.sizeout);
            y = zeros(this.sizein);
            fx = fft2(x);
            
            if this.useincident
                for iz = 1:this.Nz
                    y(:,:,iz) =   conj(this.Fincident(:,:,iz)) .* ifft2( exp( complex(0,-this.z(iz) .* this.F)) .* fx );
                end
            else
                for iz = 1:this.Nz
                    y(:,:,iz) =   ifft2( exp( complex(0,-this.z(iz) .* this.F)) .* fx );
                end
            end
        end
        function y = HtH(this,x)
            assert( isequal(size(x),this.sizein),  'x does not have the right size: [%d, %d]',this.sizein);
            y = zeros(this.sizein);
            fx = complex(zeros(this.sizeout));
            if this.useincident
                for iz = 1:this.Nz
                    fx = fx+     exp( complex(0.,this.z(iz) .* this.F)) .*  fft2(x(:,:,iz).* this.Fincident(:,:,iz));
                end
                for iz = 1:this.Nz
                    y(:,:,iz) =   conj(this.Fincident(:,:,iz)) .* ifft2( exp( complex(0,-this.z(iz) .* this.F)) .* fx );
                end
            else
                for iz = 1:this.Nz
                    fx = fx+     exp( complex(0.,this.z(iz) * this.F)) .*  fft2(x(:,:,iz));
                end
                for iz = 1:this.Nz
                    y(:,:,iz) =   ifft2( exp( complex(0,-this.z(iz) .* this.F)) .* fx );
                end
            end
        end
    end
end

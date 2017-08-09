classdef MapComposition < Map
    % MapComposition : Composition of Maps
    % $$ \\mathrm{H}(\\mathrm{x}) = \\mathrm{H}_1 \\left( \\mathrm{H}_2(\\mathrm{x}) \\right) $$
    %
    % :param H1:  left hand side :class:`Map` (or a scalar)
    % :param H2:  right hand side :class:`Map`
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2017
    %     M. McCann michael.mccann@epfl.ch
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    
    properties
        H1;           % Left hand side Map
        H2;           % Right hand side Map
        isH1Scalar=0; % True if H1 is scalar
    end
    %% Constructor
    methods
        function this = MapComposition(H1,H2)
            this.name='MapComposition';
            assert((isa(H1,'Map') || (isnumeric(H1) && isscalar(H1))),'H1 have to be a Map object or a scalar');
            assert(isa(H2,'Map'),'H2 have to be a Map');
            if ~(isnumeric(H1) && isscalar(H1))
                assert(isequal(H2.sizeout,H1.sizein),'H1 sizein is inconsistent with H2 sizeout');
                assert(~H2.isComplexOut || (H1.isComplexIn==H2.isComplexOut),'isComplexOut of H2 is inconsistent with isComplexIn of H1');
            end
            this.H1 = H1;
            this.H2 = H2;
            if isnumeric(this.H1) && isscalar(this.H1), this.isH1Scalar=1; end;
            % Set Properties
            if this.isH1Scalar    % if H1 scalar
                % norm
                if this.H2.norm ~= -1, this.norm = this.H1 * this.H2.norm; end
                % isInvertible
                if this.H1~=0 && this.H2.isInvertible, this.isInvertible=true; end;
                % isDifferentiable
                if this.H2.isDifferentiable, this.isDifferentiable=true; end
                % isComplexIn/Out
                this.isComplexIn=H2.isComplexIn;
                this.isComplexOut=H2.isComplexOut || ~isreal(this.H1);
                % sizein/out
                this.sizein=H2.sizein;
                this.sizeout=H2.sizeout;
            else  % if H1 Map
                % norm
                if this.H1.norm ~= -1 && this.H2.norm ~= -1,
                    this.norm = this.H1.norm * this.H2.norm;
                else
                    this.norm=-1;
                end
                % isInvertible
                if this.H1.isInvertible && this.H2.isInvertible, this.isInvertible=true;end
                % isDifferentiable
                if this.H1.isDifferentiable && this.H2.isDifferentiable, this.isDifferentiable=true; end
                this.isComplexIn=H2.isComplexIn;
                this.isComplexOut=H1.isComplexOut;
                % sizein/out
                this.sizein=H2.sizein;
                this.sizeout=H1.sizeout;
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - applyInverse_(this,y)
	methods (Access = protected)
		function y = apply_(this, x)
            % Reimplemented from :class:`Map`  
            if this.isH1Scalar
                y = this.H1*(this.H2.apply(x));
            else
                y = this.H1.apply(this.H2.apply(x));
            end
        end
        function x = applyJacobianT_(this, y, v)
            % Reimplemented from :class:`Map`
            if this.isDifferentiable
                if this.isH1Scalar
                    x=this.H1*this.H2.applyJacobian(v,y);
                else
                    x=this.H2.applyJacobian(v,this.H1.applyJacobian(this.H2.apply(v),y));
                end
            end
        end
        function x = applyInverse_(this, y)
            % Reimplemented from :class:`Map`
            if this.isInvertible
                if this.isH1Scalar
                    x=this.H2.applyInverse(y/this.H1);
                else
                    x=this.H2.applyInverse(this.H1.applyInverse(y));
                end
            else
                x = applyInverse_@Map(y);
            end
        end
    end	
end
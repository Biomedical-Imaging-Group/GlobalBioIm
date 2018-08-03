classdef MapMultiplication < Map
    % MapMultiplication: Multiplication of Maps
    % $$ \\mathrm{H}(\\mathrm{x}) =  \\mathrm{H}_1(\\mathrm{x}) \\times \\mathrm{H}_2(\\mathrm{x}) $$
    %
    % :param Map1:  :class:`Map` object
    % :param Map2:  :class:`Map` object
    %
    % **Example** H=MapMultiplication(Map1,Map2)
    %
    % See also :class:`Map`
    
    %%    Copyright (C) 2018
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
    
    %% Properties
    % - Public 
    properties (SetObservable, AbortSet)
        M1;
        M2;
    end
    
    %% Constructor and handlePropEvents method
    methods
        function this = MapMultiplication(Map1,Map2)
            % Listeners to PostSet events
            addlistener(this,'M1','PostSet',@this.handlePropEvents);
            addlistener(this,'M2','PostSet',@this.handlePropEvents);
            % Basic properties
            this.M1=Map1;
            this.M2=Map2;
            this.isInvertible=false;
            this.sizein = this.M1.sizein;
            this.sizeout = this.M2.sizeout;
            % Listeners to modified events (for properties which are classes) 
            addlistener(this.M1,'modified',@this.handleModifiedM1);
            addlistener(this.M2,'modified',@this.handleModifiedM2);
        end
        function handleModifiedM1(this,~,~) % Necessary for properties which are objects of the Library
            sourc.Name='M1'; handlePropEvents(this,sourc);
        end
        function handleModifiedM2(this,~,~) % Necessary for properties which are objects of the Library
            sourc.Name='M2'; handlePropEvents(this,sourc);
        end
        function handlePropEvents(this,src,~)
            % Reimplemented from parent class :class:`Map`
            switch src.Name
                case {'M1','M2'}
                    if ~isempty(this.M1) &&  ~isempty(this.M2)  % Because in the constructor they are set one after the other.
                        assert((isa(this.M1,'Map') && isa(this.M2,'Map')),'Properties M1 and M2 should be Map objects');
                        this.isDifferentiable= this.M1.isDifferentiable*this.M2.isDifferentiable;
                        assert(cmpSize(this.M1.sizein,this.M2.sizein),'Properties M1 and M2 should have consistent  sizein') ;
                        assert(cmpSize(this.M1.sizeout,this.M2.sizeout),'Properties M1 and M2 should have consistent sizeout');
                        this.name=sprintf('MapMultiplication( %s ; %s )',this.M1.name,this.M2.name);
                    end
            end
            % Call mother classes at this end (important to ensure the
            % right execution order)
            handlePropEvents@Map(this,src);
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyJacobianT_(this, y, v)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function y = apply_(this,x) 
            % Reimplemented from :class:`Map`   
            y = (this.M1*x).*(this.M2*x);
        end  
        function x = applyJacobianT_(this, y, v)
            % Reimplemented from :class:`Map`   
            
            x=(this.M1*v).*this.M2.applyJacobianT(y,v) + (this.M2*v).*this.M1.applyJacobianT(y,v);
        end     
        function M = makeComposition_(this,G)
            % Reimplemented from :class:`Map`  
            M=MapMultiplication(this.M1*G,this.M2*G);
        end
    end
end


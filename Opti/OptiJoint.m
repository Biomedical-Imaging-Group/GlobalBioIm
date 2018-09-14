classdef OptiJoint < Opti
    % Abstract class for optimization algorithms to minimize :class:`Cost` objects
    %
    %
    % See also :class:`Opti` :class:`OutputOptiJoint` :class:`Cost`
    
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
    % - Readable
    properties (SetAccess = protected,GetAccess = public)
        OptiX
        OptiZ
        cumIter % cumulative iterations
    end
    
    %% Constructor
    methods
        function this = OptiJoint(OptiX,OptiZ)
            this.OptiX=OptiX;
            this.OptiZ=OptiZ;
            this.name=['OptiJoint(X->',OptiX.name,',Z->',OptiZ.name,')'];
            this.OutOp=OutputOptiJoint(this.OptiX.OutOp,this.OptiZ.OutOp);
        end
    end
    
    %% Methods
    methods
        function run(this,x0,z0)
            % Run the algorithm.
            %
            % :param x0: initial point in \\(\\in X\\), if no argument restarts from the current value :attr:`xopt`.
            %
            % **note**: this method does not return anything, the result being stored in public attribute :attr:`xopt`.
            if(nargin==1)
                assert(~isempty(this.OptiX.xopt),'Missing starting point x0');
                assert(~isempty(this.OptiZ.xopt),'Missing starting point z0');
                x0=[];z0=[];
            end
            this.OptiX.initialize(x0);
            this.OptiZ.initialize(z0);
            this.OptiX.verbose=false;
            this.OptiZ.verbose=false;
            tstart=tic;
            this.OutOp.init();
            this.niter=0;
            this.cumIter=0;
            this.starting_verb();
            this.endingMessage= ['Maximum number of outer iterations reached: ', num2str(this.maxiter)];
            while (this.niter<this.maxiter)
                this.niter=this.niter+1;
                % - Iterations over X
                kk=0;
                while kk<this.OptiX.maxiter
                    flag=this.OptiX.doIteration();
                    if(flag == this.OPTI_NEXT_IT)
                        this.niterUp();
                        kk=kk+1;
                        % - Call OutputOpti object
                        if this.OptiX.ItUpOut>0 && (((mod(kk,this.OptiX.ItUpOut)==0)|| (kk==1)))
                            this.OptiX.OutOp.update(this.OptiX);
                            this.OutOp.update(this);
                        end
                        this.OptiX.xold=this.OptiX.xopt;
                    elseif  flag==this.OPTI_STOP
                        break;
                    end
                end
                this.upXtoCostZ();
                % - Iterations over Z
                kk=0;
                while (kk<this.OptiZ.maxiter)
                    flag=this.OptiZ.doIteration();
                    if(flag == this.OPTI_NEXT_IT)
                        this.niterUp();
                        kk=kk+1;
                        % - Call OutputOpti object
                        if this.OptiZ.ItUpOut>0 && (((mod(kk,this.OptiZ.ItUpOut)==0)|| (kk==1)))
                            this.OptiZ.OutOp.update(this.OptiZ);
                            this.OutOp.update(this);
                        end
                        this.OptiZ.xold=this.OptiZ.xopt;
                    elseif  flag==this.OPTI_STOP
                        break;
                    end
                end
                this.upZtoCostX();
%                 if ~isempty(this.OptiX.OutOp.xtrue)
%                     disp(['Iter: ',num2str(this.niter),' --- SNR coefs: ',num2str(OptSNR(this.OptiX.OutOp.xtrue,this.OptiX.xopt)),' --- Max error pos: ',num2str(max(sqrt(sum((this.OptiZ.xopt-this.OptiZ.OutOp.xtrue).^2,1))))]);
%                 else
%                     disp(['Iter: ',num2str(this.niter),' --- Max norm pos: ',num2str(max(sqrt(sum((this.OptiZ.xopt).^2,1))))]);
%                 end               
            end
            this.time=toc(tstart);
            this.ending_verb();
        end
        function niterUp(this)
            this.cumIter=this.cumIter+1;
            this.OptiX.niter=this.cumIter;
            this.OptiZ.niter=this.cumIter;
        end
        function upXtoCostZ(this)
            error('Method upZtoCostX not implemented');
        end
        function upZtoCostX(this)
            error('Method upZtoCostX not implemented');
        end
        function starting_verb(this)
            % Generic method to display a starting message in verbose mode.
            
            if this.verbose
                fprintf('---> Start %s ... \n',this.name);               
            end
        end
        function ending_verb(this)
            % Generic method to display a ending message in verbose mode.
            
            if this.verbose
                disp(this.endingMessage);
                fprintf('... Optimization finished \nElapsed time (s): %4.2d (%i outer iterations). \n',this.time, this.niter);
            end
        end
    end
end
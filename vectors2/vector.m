classdef vector<handle
    
% examples    
%
% v = vector(1,2); or v = vector([1,2]);
%
% v.plot - plot vector with no properties
% v.plot(3); line width= 3; default  line width= 0.5;
%  v.plot(3,2); plot vector  line width= 3 and color 2 from the standard Matlab color palette
%  v.plot(3,2,'-.'); %  v.plot(3,2);
%
% addPlot(v1,v2) - addition of vectors on the plane
%
% minPlot(v1, v2) - subtraction of vectors in the plane
%
% multScalarPlot(v,a) - multiplying a vector by a scalar
% 
%  linComb(v1, v2,a1,a2) - linear combination
% 
% annLinComb(v1,v2) - the animation is a linear combination of 2 vectors
    
    properties
        x
        lw
        color
        r=1
    end
    properties (Constant)
          m = [0    0.4470    0.7410;  
                                    0.8500    0.3250    0.0980;  
                                    0.9290    0.6940    0.1250; 
                                    0.4940    0.1840    0.5560; 
                                    0.4660    0.6740    0.1880; 
                                    0.3010    0.7450    0.9330; 
                                    0.6350    0.0780    0.1840];
    end
    
    methods
        function obj = vector(X1, X2)
            if nargin == 1
                if size(X1,2)~=2
                    error('the size of vectors need to be Nx2')
                end
             obj.x = X1;     
            elseif nargin == 2
                obj.x(1) = X1; obj.x(2) = X2;
            end
        end
        
        function [p, p1, p2] = plot(obj,LW,COLOR,t)
            
              if nargin ==1
                obj.lw = 2;
            else 
                obj.lw = LW;
              end
            
               if nargin >2
                    if isnumeric(COLOR)
                        obj.color = obj.m(COLOR,:);
                    else
                         obj.color = COLOR;
                    end
              else
                obj.color = 'b';
               end
              
               if nargin<4
                   t = '-';
               end
              
            ax = gca;
            grid on
            hold on
            h = false;
            if isempty(ax.Children)
                S = max(abs(obj.x(:))); 
                S = S*1.2;
                h = true;
            else 
                S = max(ax.XLim(2),max(abs(obj.x(:))));
                if ax.XLim(2)<max(abs(obj.x(:)))
                    h = true;
                end
            end
               obj.r=S/8;
            ax.XLim = [-S S]; ax.YLim = [-S S];
            if h
                plot([0 0], [-S S],'-.k')
                 plot([-S S],[0 0],'-.k')
            end
             if  obj.x(1)~=0 || obj.x(2)~=0   
            if nargout>0
                p = plot([0,obj.x(1)], [0 obj.x(2)],t,'LineWidth',obj.lw);

                   p.Color = obj.color;         
            else
                plot([0,obj.x(1)], [0 obj.x(2)],t,'LineWidth',obj.lw,'Color',obj.color);
            end

                if obj.x(1)==0 && obj.x(2)>0
                    x1 = 0; y1 = -obj.r;
                elseif  obj.x(1)==0 && obj.x(2)<0
                    x1 = 0; y1 = obj.r;
                 elseif  obj.x(1)>0 && obj.x(2)==0  
                     x1=-obj.r; y1 = 0;
                 elseif  obj.x(1)<0 && obj.x(2)==0
                    x1 = obj.r; y1 = 0;
                 elseif  obj.x(1)~=0 && obj.x(2)~=0       
                    k = obj.x(2)/obj.x(1); x1 = obj.r/sqrt(1+k^2); y1 = k*obj.r/sqrt(1+k^2);
                    if obj.x(1)*obj.x(2) >0        
                     x1 = -sign(obj.x(1))*x1; y1 = -sign(obj.x(2))*y1;
                    else
                        if sign(x1) == sign(obj.x(1))
                             x1 = -x1;
                        end
                        if sign(y1) == sign(obj.x(2))
                             y1 = -y1;
                       end
                    end
                end
                
                t = pi/180*10;
                 A1 =  [cos(t), - sin(t); sin(t), cos(t)];
                 A2 =  [cos(t),  sin(t); -sin(t), cos(t)];

                 X1 = A1*[x1;y1]; X2 = A2*[x1;y1];

                 X1 = X1 + [obj.x(1); obj.x(2)]; X2 = X2 + [obj.x(1); obj.x(2)];
                if nargout>0
                    p1 =  plot([obj.x(1),X1(1)], [obj.x(2), X1(2)],'Color',p.Color,'LineWidth',obj.lw);
                     p2 = plot([obj.x(1),X2(1)], [obj.x(2), X2(2)],'Color',p.Color,'LineWidth',obj.lw);
                else
                    plot([obj.x(1),X1(1)], [obj.x(2), X1(2)],'Color',obj.color,'LineWidth',obj.lw);
                    plot([obj.x(1),X2(1)], [obj.x(2), X2(2)],'Color',obj.color,'LineWidth',obj.lw);
                end
              else
                  plot(0,0,'ob')
             end
        end
        
        function addPlot(obj1,obj2)
            obj= vector(obj1.x + obj2.x);
            figure
            obj1.plot(2,1);
            obj2.plot(2,2);
            obj.plot(2,3);
            ax = gca;
            S = max(max(abs(obj.x(:))),ax.XLim(2));
            if  max(abs(obj.x(:))) >= ax.XLim(2)
                S = S*1.1;
                
                plot([0 0], [-S S],'-.k')
             plot([-S S],[0 0],'-.k')
            end
            ax.XLim = [-S S]; ax.YLim = [-S S];
            
            plot([obj1.x(1) obj1.x(1)+obj2.x(1)], [obj1.x(2) obj1.x(2)+obj2.x(2)],'k--')
            plot([obj2.x(1) obj1.x(1)+obj2.x(1)], [obj2.x(2) obj1.x(2)+obj2.x(2)],'k--')
            
            
        end
        
        function minPlot(obj1, obj2)
            obj= vector(obj1.x - obj2.x);
             figure
            obj1.plot(2,1)
            obj2.plot(2,2)
            obj.plot(2,3)
            minObj = vector(-obj2.x);
            minObj.plot(2,2,'--')
            
            plot([minObj.x(1) minObj.x(1)+obj1.x(1)], [minObj.x(2) minObj.x(2)+obj1.x(2)],'k--')
            plot([obj1.x(1) minObj.x(1)+obj1.x(1)], [obj1.x(2) minObj.x(2)+obj1.x(2)],'k--')
        end
        
        function multScalarPlot(obj,a)
            figure
            plot(obj,3,1)
            obj1 = vector(obj.x*a);
           obj1.plot(2,2)
        end
        
        function linComb(vec1, vec2,a1,a2)
            ax = gca; hold on
            x1 = vec1.x*a1; x2 = vec2.x*a2;
            vecT1 = vector(x1); vecT2 = vector(x2);
            S = max(abs(([x1, x2 vec1.x vec2.x])));
            if ~isempty(ax.Children)
                S = max(S,ax.XLim(2));
            end
            vecT1.plot(0.5,'k')
            vecT2.plot(0.5,'k')
           % plot(ax, [0 x1(1)], [0, x1(2)],'k')
            % plot(ax, [0 x2(1)], [0, x2(2)],'k')
            vec1.plot(3,1)
            vec2.plot(3,2)
            
            vec = vector(x1 + x2);
            plot(vec,3,4)
            S = max(S,max(abs(vec.x)));
            ax.XLim = [-S S]; ax.YLim = [-S S];
            
            plot([vecT1.x(1) vecT1.x(1)+vecT2.x(1)], [vecT1.x(2) vecT1.x(2)+vecT2.x(2)],'k--')
            plot([vecT2.x(1) vecT1.x(1)+vecT2.x(1)], [vecT2.x(2) vecT1.x(2)+vecT2.x(2)],'k--')
            
            
        end
        
        function annLinComb(vec1,vec2)
           % a1 = linspace(-2,2); a2 = linspace(-2,2);
           A = [vec1.x(1) vec2.x(1); vec1.x(2) vec2.x(2)]^-1;
           X = funcAnimate();
           A = A*X; a1 = A(1,:); a2 = A(2,:);
            ax = gca; hold on
             vec1.plot(3,1);   vec2.plot(3,2);
             S = max(abs([vec1.x(:); vec2.x(:)]));
            ax.XLim = [-S*1.5, S*1.5]; ax.YLim = [-S*1.5, S*1.5];
            
            v1 = vector(vec1.x*a1(1)); v2 = vector(vec2.x*a2(1)); v3 =vector(v1.x+v2.x);
            av1 = plot([v1.x(1) v1.x(1)+v2.x(1)], [v1.x(2) v1.x(2)+v2.x(2)],'k--');
            av2 = plot([v2.x(1) v1.x(1)+v2.x(1)], [v2.x(2) v1.x(2)+v2.x(2)],'k--');
            
            [pv11, pv12, pv13] = plot(v1,0.5,'k');  [pv21,pv22, pv23] = plot(v2,0.5,'k');
            [p31, p32,p33] = plot(v3,3,3);
            
            plot(v3.x(1),v3.x(2),'k')
            function X = funcAnimate()
                fi = linspace(0,8*pi, 400);
                R = cos(9*fi/4)+7/3;
                yy = 1.5*R.*sin(fi); xx = 1.5*R.*cos(fi);
                X = [xx;yy];
             end
            
            for i = 1:numel(a1)
                v1 = vector(vec1.x*a1(i)); v2 = vector(vec2.x*a2(i)); v3 =vector(v1.x+v2.x);
                plot(v3.x(1),v3.x(2),'k.')
                av1.XData = [v1.x(1) v1.x(1)+v2.x(1)]; av1.YData = [v1.x(2) v1.x(2)+v2.x(2)];
                av2.XData = [v2.x(1) v1.x(1)+v2.x(1)]; av2.YData = [v2.x(2) v1.x(2)+v2.x(2)];
                delete(pv11),  delete(pv12),  delete(pv13),  delete(pv21),  delete(pv22),  delete(pv23),  delete(p31),  delete(p32),  delete(p33);
                [pv11, pv12, pv13] = plot(v1,0.5,'k');  [pv21,pv22, pv23] = plot(v2,0.5,'k');
                      [p31, p32,p33] = plot(v3,3,3);
                      drawnow
                
            end
            
        end
            
    end
    

end












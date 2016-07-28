function plot_arm_sim(sd,params,varargin)

if nargin == 3 %provided tuning curve data
    % will use this to sort the raster by cell PD
    tc = varargin{1};
    use_model = 'muscle';
else
    tc = [];
end

if length(sd) > 1
    error('Only provide one trial, at least for now');
end

L1 = params.L1; % length of upper arm in m
L2 = params.L2; % length of lower arm in m
comp_blocks = params.comp_blocks;

th = sd.kin.(params.kin_model).angles;
dth = sd.kin.(params.kin_model).dangles;
ddth = sd.kin.(params.kin_model).ddangles;
p = sd.kin.(params.kin_model).pos;
v = sd.kin.(params.kin_model).vel;
T = sd.torques;
m = sd.muscles;

sd.angles = th;
sd.dangles = dth;
sd.ddangles = ddth;
sd.pos = p;
sd.vel = v;
sd.torques = T;
sd.muscles = m;

switch lower(params.type)
    case 'video'
        
        figure;
        for t = 1:size(sd.torques,1)
            hold all;
            % plot circle at endpoint
            plot(L1*cos(sd.angles(t,1)) + L2*cos(sd.angles(t,1)+sd.angles(t,2)),L1*sin(sd.angles(t,1)) + L2*sin(sd.angles(t,1)+sd.angles(t,2)),'ko','LineWidth',3);
            plot(sd.pos(t,1),sd.pos(t,2),'ro','LineWidth',2);
            
            % plot limb segments
            plot([0,L1*cos(sd.angles(t,1))],[0,L1*sin(sd.angles(t,1))],'k-','LineWidth',2);
            plot(L1*cos(sd.angles(t,1))+[0,L2*cos(sd.angles(t,1)+sd.angles(t,2))],L1*sin(sd.angles(t,1))+[0,L2*sin(sd.angles(t,1)+sd.angles(t,2))],'k-','LineWidth',2);
            set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-0.15,0.15],'YLim',[0 0.30]); axis('square');
            pause(10/1000);
            clf;
        end
        
    case 'freeze_video'
        
        % build color gradient vector
        %plot_colors = [(0:1/size(sd.angles,1):1)', zeros(size(0:1/size(sd.angles,1):1))', (1:-1/size(sd.angles,1):0)'];
        plot_colors = [(0:1/(params.idx_end-params.idx_start):1)', zeros(size(0:1/(params.idx_end-params.idx_start):1))', (1:-1/(params.idx_end-params.idx_start):0)'];
        
        figure;
        %for t = 1:params.resolution:size(sd.angles,1)
        for t = params.idx_start:params.resolution:params.idx_end
            hold all;
            % plot circle at endpoint
            plot(L1*cos(sd.angles(t,1)) + L2*cos(sd.angles(t,1)+sd.angles(t,2)), L1*sin(sd.angles(t,1)) + L2*sin(sd.angles(t,1)+sd.angles(t,2)),'o','LineWidth',3,'Color',plot_colors(t-(params.idx_start-1),:));
            plot(p(t,1),p(t,2),'o','LineWidth',2,'Color','k');
            
            % plot limb segments
            plot([0,L1*cos(sd.angles(t,1))],[0,L1*sin(sd.angles(t,1))],'k-','LineWidth',2,'Color',plot_colors(t-(params.idx_start-1),:));
            plot(L1*cos(sd.angles(t,1))+[0,L2*cos(sd.angles(t,1)+sd.angles(t,2))], L1*sin(sd.angles(t,1))+[0,L2*sin(sd.angles(t,1)+sd.angles(t,2))],'k-','LineWidth',2,'Color',plot_colors(t-(params.idx_start-1),:));
        end
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-0.21,0.21],'YLim',[-0.03 0.21]); axis('square');
        set(gca,'YTickLabel',cellfun(@(x) 100*(str2num(x) - params.origin_pos(2)/100),get(gca,'YTickLabel')))
        set(gca,'XTickLabel',100*cellfun(@(x) str2num(x),get(gca,'XTickLabel')))
        xlabel('X Position (cm)','FontSize',14);
        ylabel('Y Position (cm)','FontSize',14);
        
        
    case 'time_signals'
        
        t = params.dt*(1:size(p,1))';
        figure('Position',[200 200 400 600]);
        for i = 1:length(params.signals)
            if strcmpi(params.signals{i},'muscles') % normalize muscle activations
                subplot(length(params.signals)+1,1,i);
                plot(repmat(t,1,size(sd.(params.signals{i}),2)),sd.(params.signals{i})./repmat(params.M_max,size(sd.(params.signals{i}),1),1),'LineWidth',2);
            elseif length(params.signals{i}) > 6 && strcmpi(params.signals{i}(end-6:end),'neurons') % do raster
                subplot(length(params.signals)+1,1,i:i+1);
                if ~isempty(tc) % sort by PD and only include tuned cells
                    pds = tc(1).(use_model).tc(:,3);
                    [~,I] = sort(pds);
                    
                    is_tuned = all(prctile(tc(comp_blocks(1)).(use_model).rs,[2.5 97.5],2) > 0.1,2) & ...
                        all(prctile(tc(comp_blocks(2)).(use_model).rs,[2.5 97.5],2) > 0.1,2) & ...
                        all(prctile(tc(comp_blocks(3)).(use_model).rs,[2.5 97.5],2) > 0.1,2);% & ...
%                         angleDiff(tc(comp_blocks(1)).(use_model).cb{3}(:,1),tc(comp_blocks(1)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
%                         angleDiff(tc(comp_blocks(2)).(use_model).cb{3}(:,1),tc(comp_blocks(2)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
%                         angleDiff(tc(comp_blocks(3)).(use_model).cb{3}(:,1),tc(comp_blocks(3)).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
                    
                    I = I(is_tuned);
                else
                    I = 1:params.num_neurons;
                end
                imagesc(t-params.dt/2,1:length(I),~sd.(params.signals{i})(:,I)',[0 1]);
                colormap('gray');
            else % just plot it
                subplot(length(params.signals)+1,1,i);
                plot(repmat(t,1,size(sd.(params.signals{i}),2)),sd.(params.signals{i}),'LineWidth',2);
            end
            
            axis('tight');
            if isfield(params,[params.signals{i} '_lim'])
                set(gca,'YLim',params.([params.signals{i} '_lim']));
            end
            set(gca,'Box','off','TickDir','out','FontSize',14);
            ylabel(params.signals{i},'FontSize',14);
            
            set(gca,'XLim',[params.idx_start,params.idx_end].*params.dt);
        end
        xlabel('Time (sec)','FontSize',14);
        
    otherwise
        error('Type not recognized');
end
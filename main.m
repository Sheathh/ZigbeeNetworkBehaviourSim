clear;
close all;
timeSlots = 20000; % Simulation steps
trafficRates = [14,13,12,11,10,8,6,4,2];
trafficRates = [50:-5:15,trafficRates];
num_Coordinate = 1;
num_Router = 4 ;
num_EndDevice = 20;

% Number of Nodes
N = num_Coordinate + num_Router + num_EndDevice; 

nodeThroughput = zeros(timeSlots, N);

%% Node Definition
Node = struct('Role',       [], ...
              'Parent',     [], ...
              'Children',   [], ...
              'Queue',      [], ...
              'trafficRate',[], ...
              'EnvData',    [], ...
              'Position',   []);
% Define the Coordinator, Router and EndDevice Node
for i = 1:N
    if i == num_Coordinate
        % One Coordinator
        Node(i).Role = 'Coordinator';
        Node(i).Parent = 0;
    elseif i <= num_Coordinate + num_Router
        % All Routers are under Coordinator
        possible_parents = 1:(i-1);
        Node(i).Role = 'Router';
        Node(i).Parent = possible_parents(randi(length(possible_parents))); 
        % while(Node(i).Parent == i)
        %     Node(i).Parent = randi(num_Coordinate + num_Router); 
        % end
    else
        % Random parent node for the EndDevice
        Node(i).Role = 'EndDevice';
        Node(i).Parent = randi([num_Coordinate + 1, num_Coordinate + num_Router]); 
    end

    % Define the Children nodes, it describes the topology of the system
    Node(i).Children = [];
    % Define the packtes Queue for each node
    Node(i).Queue = [];

    % Only EndDevice node generates data packet 
    if strcmp(Node(i).Role,'EndDevice')
        % Randomly distrubute trafficrate for all the EndDevice
        % nodes.
        Node(i).trafficRate = trafficRates(ceil(randi(size(trafficRates))));
    else
        % Other nodes DO NOT generate data packet
        Node(i).trafficRate = Inf;
    end

        % Define the environment data, here choose temprature in degrees
    if strcmp(Node(i).Role,'EndDevice')
        Node(i).EnvData = 20 + 10*rand();
    else
        Node(i).EnvData = []; % For Coordinator and Router, it do not generate data
    end
    Node(i).Position = [rand()*100, rand()*100]; % Random position distribution
end
% Define the Children according to Parents
for i = num_Coordinate + 1:1:N
    index_p = Node(i).Parent;
    Node(index_p).Children = [Node(index_p).Children, i];
end


%% Packet definition
packetTemplate = struct('SrcID',        [], ...
                        'DestID',       [], ...
                        'Timestamp',    [], ...
                        'Data',         []);


% 延迟存储变量
delays_all = cell(1, length(trafficRates));
meanDelays = zeros(1, length(trafficRates));
maxDelays = zeros(1, length(trafficRates));
throughput = zeros(1, length(trafficRates));


%% Simulation for the Zigbee Behaviour
for tr = 1:length(trafficRates)
    % Set current TrafficRate
    currentTrafficRate = trafficRates(tr);
    for i = (1 + num_Router + 1):N 
        Node(i).trafficRate = currentTrafficRate;
    end
    
    delays_current = [];
    % Simulation step
    for t = 1:timeSlots
        % Data Generation
        for i = (1 + num_Router + 1):N
            if strcmp(Node(i).Role, 'EndDevice')
                % Generate data according to the traffic rate
                if rand() < (1 / Node(i).trafficRate)
                    packet = packetTemplate;
                    packet.SrcID = i;
                    packet.DestID = 1; % Destination: Coordinator
                    packet.Timestamp = t;
                    packet.Data = Node(i).EnvData;
                    Node(i).Queue = [Node(i).Queue, packet];
                end
            end
        end
        
        % 2. Data trasfer, from EndDevice to Router then to the Coordinator
        for i = N:-1:2
            if ~isempty(Node(i).Queue)
                % First packet in the Queue
                pkt = Node(i).Queue(1);
                Node(i).Queue(1) = []; % Dequeue

                nodeThroughput(t, i) = nodeThroughput(t, i) + 1;

                p = Node(i).Parent;

                % Put the packet in the end of the Parent node's queue
                Node(p).Queue = [Node(p).Queue, pkt];
            end
        end
        
        % 3. Coordinator
        if ~isempty(Node(1).Queue)
            for idx = 1:length(Node(1).Queue)
                arrival_pkt = Node(1).Queue(idx);
                delay = t - arrival_pkt.Timestamp;
                delays_current = [delays_current, delay];
            end

            numReceived = length(Node(1).Queue);
            nodeThroughput(t, 1) = nodeThroughput(t, 1) + numReceived;

            throughput(tr) = throughput(tr) + length(Node(1).Queue);
            Node(1).Queue = [];
        end
    end
    
    % Record
    if ~isempty(delays_current)
        maxDelays(tr) = max(delays_current);
        meanDelays(tr) = mean(delays_current);
    else
        maxDelays(tr) = 0;
        meanDelays(tr) = 0;
    end
    delays_all{tr} = delays_current;
    fprintf('Traffic scenario: 1 pkt per %d slots, Max Delay = %.2f, Mean Delay = %.2f, Throughput = %d\n', ...
             currentTrafficRate, maxDelays(tr), meanDelays(tr), throughput(tr));

    % Reset the Queue
    for i = 1:N
        Node(i).Queue = [];
    end
end

%% Visualization

% ========== Topology Visulization ========== %
edges = [];
for i = 2:N
    p = Node(i).Parent;
    if p > 0
        edges = [edges; p, i];
    end
end

% Digraph
G = digraph(edges(:,1), edges(:,2));

figure('Name','Tree Topology (Logical)','NumberTitle','off');
P = plot(G, ...
    'NodeColor','k', ...  
    'MarkerSize', 4, ...
    'LineWidth',1.5);
title('Zigbee Tree Network - Logical Topology');
hold on;

coordIdx  = find(strcmp({Node.Role}, 'Coordinator'));
routerIdx = find(strcmp({Node.Role}, 'Router'));
enddevIdx = find(strcmp({Node.Role}, 'EndDevice'));

highlight(P, coordIdx, ...
    'NodeColor','r', ...   
    'MarkerSize',8); 
highlight(P, routerIdx, ...
    'NodeColor','g', ...   
    'MarkerSize',6);
highlight(P, enddevIdx, ...
    'NodeColor','b', ...   
    'MarkerSize',4);

plotCoord  = plot(nan, nan, 'or', 'MarkerSize', 8,'MarkerFaceColor','r', 'DisplayName', 'Coordinator');
plotRouter = plot(nan, nan, 'og', 'MarkerSize', 6,'MarkerFaceColor','g', 'DisplayName', 'Router');
plotEnddev = plot(nan, nan, 'ob', 'MarkerSize', 4,'MarkerFaceColor','b', 'DisplayName', 'EndDevice');
legend([plotCoord, plotRouter, plotEnddev], 'Location','best');
hold off;


% ========== Node Position Visualizatoin ========== %

% Node position
positions = zeros(N, 2);
for i = 1:N
    positions(i, :) = Node(i).Position;
end

% Direction: leaf to root
H = flipedge(G);
figure('Name','Tree Topology (Physical Coordinates)','NumberTitle','off');

nodelables = 1:1:N;

H_plot = plot(H, ...
    'XData', positions(:,1), ...
    'YData', positions(:,2), ...
    'NodeLabel', nodelables, ...
    'Marker', 'none', ...    
    'LineWidth',1.5, ...
    'EdgeColor',[0 0 0]);
title('Zigbee Tree Network - Physical Positions');
hold on;

scatter(positions(coordIdx,1), positions(coordIdx,2), ...
    100, 'r', 'filled', 'DisplayName','Coordinator');
scatter(positions(routerIdx,1), positions(routerIdx,2), ...
    60, 'g', 'filled', 'DisplayName','Router');
scatter(positions(enddevIdx,1), positions(enddevIdx,2), ...
    20, 'b', 'filled', 'DisplayName','EndDevice');

legend('Location','best');
hold off;

% ========== Max delay VS Trafficrate ========== %
figure;
plot(trafficRates, maxDelays, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'r');
xlabel('Traffic Rate (1 pkt per N slots)');
ylabel('Max Delay');
title('Max Delay vs Traffic Rate');
grid on;

% ========== Histogram of delay ========== %
figure;
hold on;
colors = lines(length(trafficRates));
for tr = 1:length(trafficRates)
    histogram(delays_all{tr}, 50, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', colors(tr,:), 'LineWidth', 2);
end
hold off;
xlabel('Delay');
ylabel('Probability');
title('Delay Distribution for Different Traffic Rates');
legend(arrayfun(@(x) sprintf('1 pkt per %d slots', x), trafficRates, 'UniformOutput', false));
grid on;

% ========== Mean delay VS traffic rate ========== %
figure;
yyaxis left
plot(trafficRates, throughput, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Traffic Rate (1 pkt per N slots)');
ylabel('Throughput (Packets Received)');
title('Throughput and Mean Delay vs Traffic Rate');
grid on;

yyaxis right
plot(trafficRates, meanDelays, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'r');
ylabel('Mean Delay');
legend({'Throughput', 'Mean Delay'}, 'Location', 'best');

%% Visualization - Throughput Heatmap with Role Separation

coordEnd = num_Coordinate;
routerEnd = num_Coordinate + num_Router;

figure('Name','Node Throughput Heatmap','NumberTitle','off');
imagesc(nodeThroughput);  
colormap(jet);             
colorbar;                 
ylabel('Time Slot');
xlabel('Node Index');
title('Node Throughput (pkt/slot) Heatmap');

hold on;


line([coordEnd+0.5, coordEnd+0.5], ylim, 'Color','k','LineStyle','--','LineWidth',1.5);
line([routerEnd+0.5 routerEnd+0.5], ylim, 'Color','k','LineStyle','--','LineWidth',1.5);


text(coordEnd-3.5, max(ylim)*1.07, 'Coordinator', 'VerticalAlignment','top','Color','k','FontWeight','bold');
text(routerEnd+1, max(ylim)*0.8, 'Routers', 'VerticalAlignment','top','Color','w','FontWeight','bold');
text(N-5, max(ylim)*0.95, 'EndDevices','VerticalAlignment','top','Color','w','FontWeight','bold');

hold off;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           


set(gca, 'YDir','normal');  

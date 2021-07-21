function ImBat_MCS_Visualize_DTMC(Fnorm34,c_s_Tnorm_34,c_s_T_34,flightPaths34,sim)


    if sim == 0
        % Visualize the markov chain as a graph
        Fnorm34_first10 = Fnorm34(1:10,1:10);
        Fnorm34_first5 = Fnorm34(1:5,1:5);

        % DTMC of first 10 clusters
        mc = dtmc(Fnorm34_first10); 
        figure(); 
        h = graphplot(mc);
        c = h.EdgeColor;
        h.EdgeColor = 'k';
        mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
        h.LineWidth = mcp;
        h.MarkerSize = [size(flightPaths34.id(flightPaths34.id==1),1)/10,
                        size(flightPaths34.id(flightPaths34.id==2),1)/10,
                        size(flightPaths34.id(flightPaths34.id==3),1)/10,
                        size(flightPaths34.id(flightPaths34.id==4),1)/10,
                        size(flightPaths34.id(flightPaths34.id==5),1)/10,
                        size(flightPaths34.id(flightPaths34.id==6),1)/10,
                        size(flightPaths34.id(flightPaths34.id==7),1)/10,
                        size(flightPaths34.id(flightPaths34.id==8),1)/10,
                        size(flightPaths34.id(flightPaths34.id==9),1)/10,
                        size(flightPaths34.id(flightPaths34.id==10),1)/10];


        % DTMC of first clusters that are large
        mc = dtmc(Fnorm34_first5); 
        figure(); 
        h = graphplot(mc);
        c = h.EdgeColor;
        h.EdgeColor = 'k';
        mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
        h.LineWidth = mcp;
        h.MarkerSize = [size(flightPaths34.id(flightPaths34.id==1),1)/10,
                        size(flightPaths34.id(flightPaths34.id==2),1)/10,
                        size(flightPaths34.id(flightPaths34.id==3),1)/10,
                        size(flightPaths34.id(flightPaths34.id==4),1)/10,
                        size(flightPaths34.id(flightPaths34.id==5),1)/10];


        %% Network of stereoyped to non

        mc = dtmc(c_s_Tnorm_34); 
        figure(); 
        h = graphplot(mc);
        c = h.EdgeColor;
        h.EdgeColor = 'k';
        mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*8;
        h.LineWidth = mcp;
        h.MarkerSize = [(c_s_T_34(1,1)+c_s_T_34(1,2))/10,
                        (c_s_T_34(2,1)+c_s_T_34(2,2))/10];

        disp(c_s_T_34(1,1)+c_s_T_34(1,2));
        disp(c_s_T_34(2,1)+c_s_T_34(2,2));
    else 
        % Visualize the markov chain as a graph
        Fnorm34_first10 = Fnorm34(1:10,1:10);
        Fnorm34_first5 = Fnorm34(1:5,1:5);

        % DTMC of first 10 clusters
        mc = dtmc(Fnorm34_first10); 
        figure(); 
        h = graphplot(mc);
        c = h.EdgeColor;
        h.EdgeColor = 'k';
        mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
        h.LineWidth = mcp;
        
        % DTMC of first clusters that are large
        mc = dtmc(Fnorm34_first5); 
        figure(); 
        h = graphplot(mc);
        c = h.EdgeColor;
        h.EdgeColor = 'k';
        mcp = nonzeros(reshape(mc.P,[size(mc.P,1)*size(mc.P,1),1])).*5;
        h.LineWidth = mcp;

end
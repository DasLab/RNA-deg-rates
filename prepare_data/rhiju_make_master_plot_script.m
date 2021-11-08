%% read in the data
data_files = {...
    'RollYourOwnStructure1_CloudLabMAPseq_EtOhPrecip_degradation_rates_Lib062020.csv',...
    'RollYourOwnStructure2_CloudLabMAPseq_EtOhPrecip_degradation_rates_Lib090720.csv',...
    'RT-CE_062620_P4P62HP_stability.csv',...
    'RollYourOwnStructure1_110620_Followup_ModifiedRNA_CE_degradation_rates.csv',...
    'RTPCR_degradation_rates_081120_233x.csv',...
    '081020_mRNA-deg-gel_LENGTHCUT.csv',...
    'Bioanalyzer_6x_mRNAs_101520.csv',...
    'Bioanalyzer_6x_mRNAs_PSU_5mC_102120_Revisit_MATLAB.csv',...
    'Bioanalyzer_18xcombos_varyCDS_varyUTR_mRNAs_113020.csv',...
    'Bioanalyzer_averaged_k_deg_Final24_MATLAB_2020only.csv',...
    'Bioanalyzer_2021_Stanford_S2P_unmodified.csv',...
    'Bioanalyzer_2021_Stanford_S2P_m1PSU.csv',...
    'Bioanalyzer_2021_Stanford_S2P_unmodified.csv',...
    'Bioanalyzer_2021_Stanford_S2P_m1PSU.csv',...
}
%    'RollYourOwnStructure1_081720_Followup_CE_degradation_rates.csv',...
%    'RTPCR_degradation_rates_082520_233x.csv',...
%    '100720_mRNA-deg-gel.csv',...

start_pos = []; end_pos = [];
k_deg = []; k_deg_err = [];
data_set_number = [];
RNA_type = {};
RNA_sequences_data = {};
for i = 1:length( data_files )
    x = readtable(['data/',data_files{i}],'ReadVariableNames',1 ,'Delimiter',',','VariableNamingRule','preserve');
    dx = table2cell( x );
    RNA_type = [RNA_type, dx(:,4)'];
    start_pos = [start_pos, cell2mat(dx(:,5))'];
    end_pos   = [end_pos, cell2mat(dx(:,6))'];
    k_deg = [k_deg, cell2mat(dx(:,7))'];
    k_deg_err = [k_deg_err, cell2mat(dx(:,8))'];
    RNA_sequences_data = [RNA_sequences_data, strrep(upper(dx(:,9)),'T','U')'];
    data_set_number = [data_set_number, i * ones(1,size(dx,1))];
end
L = end_pos - start_pos + 1;

    
% plot the data -- just lengths
set(figure(1),'pos',[36   553   800   800]); clf;

unique_RNA_types = unique(RNA_type,'stable');
%unique_RNA_types = setdiff(unique(RNA_type,'stable'),{'COV2 Eterna'});
for i = 1:length(unique_RNA_types)
    gp = find(strcmp( RNA_type, unique_RNA_types{i} ) );
    %plot( L(gp), k_deg(gp),'x','linew',2 ); hold on
    L_jitter = L(gp)+5*randn(1,length(gp));
    h = errorbar( L_jitter, k_deg(gp),k_deg_err(gp),'.'); hold on
    if length(gp) < 25;
        set(h,'marker','o','markerfacecolor',get(h,'markeredgecolor'));
    else
        set(h,'marker','o','markerfacecolor','w');
    end
    set(h,'color',get_color_for_master_plot(unique_RNA_types{i}))
    set(h,'marker',get_marker_for_master_plot(unique_RNA_types{i}))
end

k_cleave_LiBreaker = (7.91E-05) * 60;
correct_LiBreaker = 1/3;
x= [1:5000];
%plot(x, k_cleave_LiBreaker*x, 'k-','color',[0.9 0.9 0.9]);
plot(x, correct_LiBreaker *    k_cleave_LiBreaker*x, 'k-','color',[0.7 0.7 0.7]); hold on
plot(x, correct_LiBreaker *0.4*k_cleave_LiBreaker*x,'-','color',[0.5 0.5 0.5],'linew',1.5);
plot(x, correct_LiBreaker *0.2*k_cleave_LiBreaker*x, '-','color','k','linew',1.5 );

legend( [strrep(unique_RNA_types,'_',' '),{'1/3x Li,Breaker [AUP=1]','1/3x Li,Breaker [AUP=0.4]','1/3x Li,Breaker [AUP=0.2]'}],'location','northwest');
set(gcf, 'PaperPositionMode','auto','color','white');
xlabel( 'Length (nucleotides)');
ylabel( 'k_{deg} (1/hr)');
set(gca,'fontsize',11,'fontweight','bold','linew',1.5);


set(figure(2),'pos',[136   753   800   800]); clf;

% csv files of predictions, and RNA sequences.
prediction_file_sets = {...
    {'predictions/233x/P_UNP_233x_EternaFold.txt','predictions/233x/RNASEQUENCES_233x.txt'},...
    {'predictions/RYOS1/P_UNP_eternafold.txt','predictions/RYOS1/RNA_sequences.txt'},...
    {'predictions/RYOS2/P_UNP_EternaFold_RYOS_r2.txt','predictions/RYOS2/RNA_sequences_RYOS_R2.txt'},...
    {'predictions/P4-P6/P_UNP_EternaFold_P4P6_sequence.txt','predictions/P4-P6/P4P6_sequence.txt'},...
    {'predictions/24FinalNLuc/P_UNP_24xNLuc_EternaFold.txt','predictions/24FinalNLuc/RNASEQUENCES_24xNLuc.txt'},...
    {'predictions/18CombosNLuc/P_UNP_18CombosNLuc_EternaFold.txt','predictions/18CombosNLuc/RNASEQUENCES_18CombosNLuc.txt'},...
    {'predictions/S2P/P_UNP_S2P_EternaFold.txt','predictions/S2P/RNASEQUENCES_S2P.txt'},...
    };

RNA_sequences_pred = {};
pred          = [];
for i = 1:length( prediction_file_sets )
    prediction_file_set = prediction_file_sets{i};

    new_pred = csvread(prediction_file_set{1})';
    
    count = size(pred,2);
    pred( (size(pred,1)+1):size(new_pred,1),1:count) = 0; % reshape
    pred = [pred, zeros(size(pred,1),size(new_pred,2))]; % add some zeros for new entries
    pred(  1:size(new_pred,1), count+[1:size(new_pred,2)]) = new_pred;
    
    new_RNA_sequences = table2cell(readtable(prediction_file_set{2},'ReadVariableNames',0,'ReadRowNames',0));
    assert( length( new_RNA_sequences)  == size( new_pred, 2 ) );

    new_RNA_sequences = strrep( upper(new_RNA_sequences),'T','U');
    RNA_sequences_pred = [ RNA_sequences_pred, new_RNA_sequences'];
end
   
assert( length( RNA_sequences_pred)  == size( pred, 2 ) );

% (WOULD BE BETTER TO TRANSPOSE PRED TO MATCH PRIOR CONVENTION!!)
figure(2)
k_pred = [];
clf;
for i = 1:length(k_deg)
    k_pred(i) = NaN;
    idx = find(strcmp(RNA_sequences_pred, RNA_sequences_data{i}));
    if isempty( idx ); continue;end;
    
    pred_profile = pred(:,idx);

    U_pos = strfind( RNA_sequences_data{i},'U' );
    if ~isempty(strfind(RNA_type{i},'m1PSU'))
        U_pos = strfind( RNA_sequences_data{i},'U' );
        pred_profile( U_pos )  = 0;
        %pred_profile( U_pos-1 )= 0;
        %pred_profile( U_pos-1 )= 0.5 * pred_profile( U_pos-1 );
    elseif ~isempty(strfind(RNA_type{i},'PSU'))
        U_pos = strfind( RNA_sequences_data{i},'U' );
        pred_profile( U_pos )  = 0;
        %pred_profile( U_pos-1 )= 0;
        %pred_profile( U_pos-1 )= 0.5 * pred_profile( U_pos-1 );
    end
    k_pred(i) = sum( pred_profile([start_pos(i):end_pos(i)])' );
end

clf;
k_cleave_LiBreaker = (7.91E-05) * 60;
correct_LiBreaker = 1/3;
x= [1:5000];
for i = 1:length(unique_RNA_types)
    gp = find(strcmp( RNA_type, unique_RNA_types{i} ) );
    if length(gp) > 1000; 
        h = plot( k_pred(gp), k_deg(gp),'.','markersize',0.5); hold on
    else
        h = errorbar( k_pred(gp), k_deg(gp),k_deg_err(gp),'.'); hold on
        if length(gp) < 25;
            set(h,'marker','o','markerfacecolor',get(h,'markeredgecolor'));
        else
            set(h,'marker','o','markerfacecolor','w');
        end        
    end
    set(h,'color',get_color_for_master_plot(unique_RNA_types{i}))
    set(h,'marker',get_marker_for_master_plot(unique_RNA_types{i}))
end
plot(x, correct_LiBreaker *    k_cleave_LiBreaker*x, 'k-'); hold on
legend( [strrep(unique_RNA_types,'_',' '),{'1/3x Li,Breaker'}],'location','northwest');
%legend( [strrep(data_files,'_',' '),{'1/2x Li,Breaker'}],'location','northwest');
set(gcf, 'PaperPositionMode','auto','color','white');
xlabel( 'Summed unpaired probability, SUP (EternaFold)');
ylabel( 'k_{deg} (1/hr)');
xlim([0 2200])
ylim([-0.5 8])
set(gca,'fontsize',11,'fontweight','bold','linew',1.5);


%% Output data to final, easy-to-use .csv files
outfile = 'RNA_deg_rates.csv';
fid = fopen( outfile, 'w' );
fprintf(fid, 'data_file,RNA_type,start_pos,end_pos,k_deg,k_deg_err,k_pred_eternafold,RNA_sequence\n');
count = 0;
for i = [3:length(unique_RNA_types),1,2]
    gp = find(strcmp( RNA_type, unique_RNA_types{i} ) );
    for n = gp
        if isnan( k_pred(n) ) continue; end;
        fprintf(fid, '%s,%s,%d,%d,%8.6f,%8.6f,%8.6f,%s\n',...
            data_files{data_set_number(n)},strrep(RNA_type{n},',',' '),start_pos(n),end_pos(n),k_deg(n),k_deg_err(n),...
            correct_LiBreaker *  k_cleave_LiBreaker*k_pred(n),...
            RNA_sequences_data{n});
        count = count + 1;
    end
end
fprintf( 'Outputted %d entries to %s.\n',count,outfile );
fclose(fid);


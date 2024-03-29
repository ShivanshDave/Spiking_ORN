%% MATLAB B2021a
% ## ORN operation
% ### Effect of Concentration ( R99 F2-A )
script_fig_txn_compare_conc;
 
%% ### Effect of duration
script_fig_txn_compare_duration; 
%% ### Adaptation
script_fig_txn_compare_adaptation;

%% ## ML spikes {#ML}
% ### CaFR modulation, nK
script_fig_ML_spikes;

%% Adding ML spikes into ORN Txn
script_fig_ML_spikes_with_ORN;

%% ## Spike identification
script_fig_Spike_ID;

% Different than adding
% exportgraphics(f,'.\Report\figs\supp_fig_ML_spikes_OVERLAP.png','Resolution',300)
% without CaFR
% exportgraphics(f,'.\Report\figs\supp_fig_ML_ID_w_CaFR_1.png','Resolution',300)
%% # Results

% % ## Spiking in ORN

%% ### Response to concentration changes
 script_fig_spk_compare_conc;

%% ### Response to duration changes
script_fig_spk_compare_duration;

%% ### Response to Adaptation   << R99.F5 >>
script_fig_spk_compare_adaptation; 

%% All together
script_fig_spk_all_comps;

%% ## Optimal sniffing frequency

% % ### steady-state FR vs sniffing frequency
script_fig_sniff_freq_FRss; 


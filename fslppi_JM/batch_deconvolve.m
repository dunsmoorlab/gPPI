all_sub_args = {1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125};
day1 = {'baseline','acquisition','extinction'};
day2 = {'memory_run-01','memory_run-02','memory_run-03'};
phases = {'baseline','acquisition','extinction','memory_run-01','memory_run-02','memory_run-03'};
%rois = {'sgACC','rACC','lh_hpc','rh_hpc','lh_amyg','rh_amyg'};
rois = {'rh_hpc'}

for subi = 1:length(all_sub_args)
    fsub = sprintf('sub-FC%03d',all_sub_args{subi});
    
    for phasei = 1:length(phases)
        phase = phases{phasei};
        featdir = sprintf('/scratch/05426/ach3377/fc-bids/derivatives/model/%s/%s/%s_%s_gPPI.feat',fsub,phase,fsub,phase);
        
        for roii = 1:length(rois)
            roi = rois{roii};
            out_dir = sprintf('/scratch/05426/ach3377/fc-bids/derivatives/model/%s/%s/%s/',fsub,phase,roi);
            voi = strcat(out_dir,sprintf('%s_bold_signal.txt',roi));
            
            if any(strcmp(day2,phase))
                con_n = 8;
                con_len = 8;
            elseif strcmp(phase,'acquisition')
                con_n = 2;
                con_len = 3;
            else
                con_n = 2;
                con_len = 2;
            end
            
            for coni = 1:con_n
                con = zeros(1,con_len);
                con(coni) = 1;
                disp(con)
                
                [PPI JUNK]=fsl_ppi(featdir, voi, con, 0);
                ppi_out=PPI.ppi;
                save(strcat(out_dir,sprintf('%s_ppi.txt',PPI.name{:})), 'ppi_out', '-ascii');
            
            end

        end

    end

end
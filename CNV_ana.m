Preprocessing   = 1;
figure_vis      = 1;

% This script infers Copy Number variation (CNV) from scRNA-seq data
% inspiret by Tirosh et al. 2014. The output is a figure showing the CNV 
% profile for each cells (every row is a cell) 
% 
% dependencies: MAGIC as an imputution method
% https://github.com/KrishnaswamyLab/MAGIC
%
%                           INPUT ARGUMENTS
%                           ---------------
%   A_fil - a filtered expression matrix (cells as columns, genes as rows)
%   clustLabel - a numeric vector indicating the cluster for each cell in A
%   (cancer, immune, stroma...)
%   clusNames - a string vector indicating the cluster names for each of the
%   values of clustLabels
%   gene_names - all of the gene names of A_fil
%   gene_loc - a matrix in the length of gene_names containg 3 columns:
%   start location, end location and chr number
%   gene_name_ord - all of the gene names of gene loc
%   sample_ID_fil - sample id (for multiple batches)
%   set of genes to remove - if wanted, you can remove specific genes from
%   the analysis (for example: differentially expressed genes)
%
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The results of this script is infered CNV from the scRNA-seq plotted 
%   on a graph

% loading data and color blind color map
try
    gene_loc;
catch
    load CNV_stuff.mat gene_loc gene_names_ord;
    load ../ZF4b/Analysis/general.mat
    load ../ZF4b/Analysis/cmap_color_blind2.mat;
    set(0,'DefaultFigureColormap',(cmap_color_blind2(1:4:end,:)));
    load ../ZF4a/working_data_ZF4a.mat
    load ../ZF4a/cluster_ZF4a.mat
    load ../ZF4b/Analysis/clone_genes.mat
end

%parameters set by user:
num_of_chr = 25; %num of chr in the genome. for zf, set as 25.
immune_clus = 2; cancer_clus = 1; % set clusters according to clustLabel
window_p = 30; % num of genes for moving average
num_of_cells = 50; % minimum number of cells expressing each gene
num_clusters = 2; % number of CNV clusters
%thersholds of CNV plot to smooth the profile
low_pass = -3;     % anything lower than this will be set as low pass   
high_pass = 3;     % anything higher than this will be set as low pass   
middle_pass_low = -1.75;  % anything higher than this but above zero will be set as 0
middle_pass_high = 1.75;  % anything lower than this but above zero will be set as 0

if (Preprocessing)
    %PREPROCESSING
    % for each gene set the average location
    % order the genes, per matrix
    % Add label for chromosome
    %make matrix gene_loc2 where columns are:
    %1. index (out of the all genes)
    %2. chromosome number
    
    %for each gene in the unorder gene_names find the index in the order
    %format
    clear gene_loc2;
    a = gene_loc(find(gene_loc(:,3)<=num_of_chr),3);
    [i,j,k] = intersect(gene_names_ord(1:length(a)),  gene_names(:,1));
    %gene index
    gene_loc2(j,1) = k;
    %chrom num
    gene_loc2(:,2) = a;
    
    % prepare the matrix - normilzation, transformation and smoothing
    try
        A_tpm;
    catch
        % normilzation
        A_tpm = median(sum(A_fil))*bsxfun(@rdivide,A_fil,sum(A_fil));
        % transformation
        A_anscombe = sqrt(A_tpm)+sqrt(A_tpm+1);  
        %smoothing
        try
            A_smooth;
        catch
            A_smooth_tpm = run_magic(A_fil',1)';
            A_smooth_trans = sqrt(A_smooth_tpm)+sqrt(A_smooth_tpm+1);
        end
        
        A_ref  = A_smooth_trans(:,find(clustLabel==immune_clus)); %these are immune cells
        A_test  = A_smooth_trans(:,find(clustLabel==cancer_clus)); %these are cancer cells
        % for each gene: shift so that reference values are zero
        A_ref_shifted = A_ref - (repmat(mean(A_ref')',1,size(A_ref,2)));
        A_test_shifted = A_test - (repmat(mean(A_ref')',1,size(A_test,2)));
        
        %use the chromosome order
        A_ref_shifted_location = A_ref_shifted (gene_loc2(:,1),:);
        A_test_shifted_location = A_test_shifted (gene_loc2(:,1),:);
        disp ('done with pre processing')
    end
    
    
    % for each cell moving average for each chromosome
    %removing genes with zero exp before moving avg
    idx = A_fil(gene_loc2(:,1),:)~=0;
    c = sum(idx,2);
    genes_with_exp = c>num_of_cells; %at least x cells
    %remove also archetype genes
    set1_genes = ismember(gene_names(gene_loc2(:,1)),clone1_genes);
    set2_genes = ismember(gene_names(gene_loc2(:,1)),clone2_genes);
    set3_genes = ismember(gene_names(gene_loc2(:,1)),clone3_genes);
    set4_genes = ismember(gene_names(gene_loc2(:,1)),clone4_genes);
    genes_to_remove = set1_genes + set2_genes + set3_genes + set4_genes;

    genes_to_keep = genes_with_exp-genes_to_remove;
  
    A_ref_shifted_location_no_zero = A_ref_shifted_location(find(genes_to_keep),:);
    A_test_shifted_location_no_zero = A_test_shifted_location(find(genes_to_keep),:);
    gene_loc2_no_zero = gene_loc2(find(genes_to_keep),:);
     
    %creating moving average for each chrom
    clear A_ref_moving A_test_moving Chrom_locations win_sizes;
    start_ma = 1;
    for chrom = 1:25
        i = find(gene_loc2_no_zero(:,2)==chrom);
        window_size = floor(length(i)/window_p);
        win_sizes(chrom)   = window_size;
        a = moving_average(A_ref_shifted_location_no_zero(i,:),window_size);
        b = moving_average(A_test_shifted_location_no_zero(i,:),window_size);
        end_ma = start_ma+window_size-1;
        A_ref_moving(:,start_ma:end_ma) = a;
        A_test_moving(:,start_ma:end_ma) = b;
        Chrom_locations(start_ma:end_ma) = chrom;
        start_ma = end_ma + 1;
    end
    disp ('done with moving average')
    line_locations = find(Chrom_locations(1:end-1)-Chrom_locations(2:end));
end

if (figure_vis)     
    A_ref_moving_zscore = zscore(A_ref_moving')';
    A_ref_moving_zscore(find(A_ref_moving_zscore<low_pass)) = low_pass;
    A_ref_moving_zscore(find(A_ref_moving_zscore>high_pass)) = high_pass;
    i = intersect(find(A_ref_moving_zscore>middle_pass_low),find(A_ref_moving_zscore<middle_pass_high)); 
    A_ref_moving_zscore(i) = 0;

    A_test_moving_zscore = zscore(A_test_moving')'; 
    A_test_moving_zscore(find(A_test_moving_zscore<low_pass)) = low_pass; 
    A_test_moving_zscore(find(A_test_moving_zscore>high_pass)) = high_pass;
    i = intersect(find(A_test_moving_zscore>middle_pass_low),find(A_test_moving_zscore<middle_pass_high)); 
    A_test_moving_zscore(i) = 0;
    
    A_ref_to_plot  = A_ref_moving_zscore;
    A_test_to_plot = A_test_moving_zscore;

    CNV_clust = cluster(linkage(pdist(A_test_to_plot),'ward'),'maxclust',num_clusters);
    %[CNV_clust,~,sumd] =  kmeans(A_test_to_plot,num_clusters,'Distance','sqeuclidean');
    [~,kk] = sort(CNV_clust);
    disp ('done with clustering')

    figure;   
    imagesc(A_test_moving_zscore(kk,:)); colorbar; 
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    for i = 1:24
        line([line_locations(i)+0.5 line_locations(i)+0.5], [0 size(A_test_moving_zscore,1)],'color','k'); hold on;
    end
    colorbar('Position',...
        [0.910763888888888 0.12035010940919 0.0111111111111118 0.155361050328228]);
    cnv_colors = brewermap(max(CNV_clust),'Set2');
    figure; imagesc(CNV_clust(kk,:)); colormap(cnv_colors)
    set(gca,'xtick',[]); set(gca,'ytick',[]);

end

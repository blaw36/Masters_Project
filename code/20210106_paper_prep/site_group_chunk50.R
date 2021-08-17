## Script to split sites into groups of arbitrary size

split_site_groups <- function(txt_file
	, grp_path = "/home/bklaw/paper_data_clean/site_splits/"
	, output_path
	, chunk_size = 50){
	site_list = read.table(paste0(grp_path, txt_file), header = T)
	num_sites = nrow(site_list)

	num_chunks = ceiling(num_sites/chunk_size)
	for(i in 1:num_chunks){
		start = ((i-1)*50) + 1
		end = i*50

		# print(start)
		# print(end)
		# print(paste0(grp_path,output_path,"chunk",i,".txt"))
		write.table(site_list[start:end,],paste0(grp_path,output_path,"chunk",i,".txt"),row.names = F)
	}
}

split_site_groups(txt_file = "sites_grp1.txt"
	, grp_path = "/home/bklaw/paper_data_clean/site_splits/"
	, output_path = "grp1/"
	, chunk_size = 50)
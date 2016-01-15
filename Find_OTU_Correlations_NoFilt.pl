# For a given diversity profile, find the corresponding epiphytes and endophytes of the "check" group and find those that correlate to the diversity

# First run export PERL5LIB=$PERL5LIB:/Users/agler/packages/Statistics-R-0.33/lib/:/Users/agler/packages/IPC-Run-0.94/lib/:/Users/agler/packages/Regexp-Common-2013031301/lib/
# export R_LIBS="/Users/agler/R-packages/"

# Filters for OTUs are:
# smr: ratio of sum of abundance to max abundance > 3.  Essentially how abundant or widespread the points are
# otusum: sum of abundance of the OTU > 50
# numpos: number of samples where the OTU shows up > 10.

use Statistics::R;

my $otutablesfile = shift; # This will be a list where each line has an OTU table followed by its mapfile, followed by its subsampling depth, and then the same for an OTU table to measure correlations against

my $otudir = "./OTU_tables/";
my $mapdir = "./Mapfiles/";
my $outdir = "./Correlation_Resubmission_Test_NoFilt/";
my $outsumdir = "./Correlation_Resubmission_Test_NoFilt/SummaryDirs/";

if(-e $outdir){}else{
    my $command1 = "mkdir $outdir";
    my $cmd1 = system('bash', '-c', "$command1") == 0
    or die "system failed: $?";
}

if(-e $outsumdir){}else{
    my $command1 = "mkdir $outsumdir";
    my $cmd1 = system('bash', '-c', "$command1") == 0
    or die "system failed: $?";
}

open TABLELIST, $otutablesfile;

while(<TABLELIST>){
    
    my $currentline = $_;
    chomp $currentline;
    
    my @splitline = split(/\t/,$currentline);
    
    my $linelen = @splitline;
    
    my $otutable1 = $splitline[0];
    my $otumapfile1 = $mapdir.$splitline[1];
    my $otumapbase1 = $splitline[1];
    my $depth1 = $splitline[2];
    
    my $otutable2 = $splitline[3];
    my $otumapfile2 = $mapdir.$splitline[4];
    my $otumapbase2 = $splitline[4];
    my $depth2 = $splitline[5];
    
    my @splitotutable1 = split(/\./, $otutable1);
    my $realbase1 = $splitotutable1[0];
    my $otutablebasename1 = $otudir.$realbase1;
    
    my @splitotutable2 = split(/\./, $otutable2);
    my $realbase2 = $splitotutable2[0];
    my $otutablebasename2 = $otudir.$realbase2;

    my @splitmapfile1 = split(/\./, $otumapbase1);
    my $realmapbase1 = $splitmapfile1[0];
    my $mapbasename1 = $mapdir.$splitmapfile1[0];
    
    my @splitmapfile2 = split(/\./, $otumapbase2);
    my $realmapbase2 = $splitmapfile2[0];
    my $mapbasename2 = $mapdir.$splitmapfile2[0];
    
    ####### Perform rarefaction and taxonomy summaries, if necessary
        
    my $rarefiedtable1 = $outsumdir.$realbase1."_".$depth1.".biom";
    my $intable1 = $otudir.$otutable1;
    
    if(-e $rarefiedtable1){}else{
    my $command1 = "single_rarefaction.py -i $intable1 -o $rarefiedtable1 -d $depth1";
    my $cmd1 = system('bash','-c',". /macqiime/configs/bash_profile.txt && $command1") == 0
    or die "system failed: $?";
    sleep 1 until -e $rarefiedtable1;   
    }
    
    my $rarefiedtable2 = $outsumdir.$realbase2."_".$depth2.".biom";
    my $intable2 = $otudir.$otutable2;
    
    if(-e $rarefiedtable2){}else{
    my $command1 = "single_rarefaction.py -i $intable2 -o $rarefiedtable2 -d $depth2";
    my $cmd1 = system('bash','-c',". /macqiime/configs/bash_profile.txt && $command1") == 0
    or die "system failed: $?";
    sleep 1 until -e $rarefiedtable2;   
    }
    
    my $summarydir1 = $outsumdir.$realbase1."_".$depth1."_summtax/";
    my $lookfor1 = $summarydir1.$realbase1."_".$depth1."_L7.txt";
    
    if(-e $summarydir1){}else{
    my $command2 = "summarize_taxa.py -i $rarefiedtable1 -o $summarydir1 -a -L 4,5,6,7";
    my $cmd2 = system('bash', '-c', ". /macqiime/configs/bash_profile.txt && $command2") == 0
    or die "system failed: $?";
    sleep 1 until -e $lookfor1;
    }

    my $summarydir2 = $outsumdir.$realbase2."_".$depth2."_summtax/";
    my $lookfor2 = $summarydir2.$realbase2."_".$depth2."_L7.txt";
    
    if(-e $summarydir2){}else{
    my $command2 = "summarize_taxa.py -i $rarefiedtable2 -o $summarydir2 -a -L 4,5,6,7";
    my $cmd2 = system('bash', '-c', ". /macqiime/configs/bash_profile.txt && $command2") == 0
    or die "system failed: $?";
    sleep 1 until -e $lookfor2;
    }

    
    my $R = Statistics::R->new();
    
    # Initialize the final output tables
    $R-> send("edgetable_4 <- t(data.frame(row.names= c(\"OrgID\",\"targetID\",\"tub_support\",\"Obs_In\",\"slope\",\"r_sq\",\"adj_r_sq\",\"p_val\",\"fdr_p_value\",\"sumtomax\",\"stm_rsqadj\", \"Sub_Support\")))");
    $R-> send("edgetable_5 <- t(data.frame(row.names= c(\"OrgID\",\"targetID\",\"tub_support\",\"Obs_In\",\"slope\",\"r_sq\",\"adj_r_sq\",\"p_val\",\"fdr_p_value\",\"sumtomax\",\"stm_rsqadj\", \"Sub_Support\")))");
    $R-> send("edgetable_6 <- t(data.frame(row.names= c(\"OrgID\",\"targetID\",\"tub_support\",\"Obs_In\",\"slope\",\"r_sq\",\"adj_r_sq\",\"p_val\",\"fdr_p_value\",\"sumtomax\",\"stm_rsqadj\", \"Sub_Support\")))");
    $R-> send("edgetable_7 <- t(data.frame(row.names= c(\"OrgID\",\"targetID\",\"tub_support\",\"Obs_In\",\"slope\",\"r_sq\",\"adj_r_sq\",\"p_val\",\"fdr_p_value\",\"sumtomax\",\"stm_rsqadj\", \"Sub_Support\")))");
    
    $R-> send("nodetable_4 <- t(data.frame(row.names= c(\"ID\",\"type\",\"EndoEpi\",\"FromWhich\",\"kingdom\",\"phyla\",\"class\",\"order\",\"family\",\"genus\",\"sum_abund\",\"max_abund\",\"num_obs\")))");
    $R-> send("nodetable_5 <- t(data.frame(row.names= c(\"ID\",\"type\",\"EndoEpi\",\"FromWhich\",\"kingdom\",\"phyla\",\"class\",\"order\",\"family\",\"genus\",\"sum_abund\",\"max_abund\",\"num_obs\")))");
    $R-> send("nodetable_6 <- t(data.frame(row.names= c(\"ID\",\"type\",\"EndoEpi\",\"FromWhich\",\"kingdom\",\"phyla\",\"class\",\"order\",\"family\",\"genus\",\"sum_abund\",\"max_abund\",\"num_obs\")))");
    $R-> send("nodetable_7 <- t(data.frame(row.names= c(\"ID\",\"type\",\"EndoEpi\",\"FromWhich\",\"kingdom\",\"phyla\",\"class\",\"order\",\"family\",\"genus\",\"sum_abund\",\"max_abund\",\"num_obs\")))");
    
    # Go ahead and load the OTU tables and mapfiles into R
    
    # Load the mapfiles for the OTU tables
    $R-> send("mapfile1_otutable1 <- read.table(\"$otumapfile1\",header=T, row.names=1, comment.char=\"\", sep=\"\t\")");
    $R-> send("mapfile2_otutable2 <- read.table(\"$otumapfile2\",header=T, row.names=1, comment.char=\"\", sep=\"\t\")");
    
    # Load the OTU tables
    my $tablegroup1;
    my $tablegroup2;
    for(my $level=4;$level<8;$level++){
        my $currenttable1 = $summarydir1.$realbase1."_".$depth1."_L".$level.".txt";
        my @splittablename = split(/_/, $realbase1);
        $tablegroup1 = $splittablename[0];
        $endoepi1 = $splittablename[3];
        $groupname1 = $splittablename[0].$splittablename[3];
        
        $R-> send("otutable1_L$level <- read.table(\"$currenttable1\", row.names=1, comment.char=\"\", sep=\"\t\", header=T)");
        $R-> send("otutable1_L$level <- t(otutable1_L$level)");
        $R-> send("log_otutable1_L$level <- log10(otutable1_L$level + 1)");
        
        my $currenttable2 = $summarydir2.$realbase2."_".$depth2."_L".$level.".txt";
        my @splittablename = split(/_/, $realbase2);
        $tablegroup2 = $splittablename[0];
        $endoepi2 = $splittablename[3];
        $groupname2 = $splittablename[0].$splittablename[3];
        
        $R-> send("otutable2_L$level <- read.table(\"$currenttable2\", row.names=1, comment.char=\"\", sep=\"\t\", header=T)");
        $R-> send("otutable2_L$level <- t(otutable2_L$level)");
        $R-> send("log_otutable2_L$level <- log10(otutable2_L$level + 1)");
        
    }
    
    my $resultsdir = $outdir.$groupname1."_vs_".$groupname2."/";
    if(-e $resultsdir){}else{
    my $command1 = "mkdir $resultsdir";
    my $cmd1 = system('bash', '-c', "$command1") == 0
    or die "system failed: $?";
    }
    
    
    # Trim the otutable mapfile to the same rows as the OTU tables (all OTU tables should be the same as differences only arose during rarefaction)
    # otu tables should only be subset of the mapfiles, so just trim the mapfile to one otu table (doesn't matter which one)
    
    $R-> send("reorder_mapfile1_otutable1 <- mapfile1_otutable1[match(row.names(log_otutable1_L4),row.names(mapfile1_otutable1)),]");
    $R-> send("nona_reorder_mapfile1_otutable1 <- reorder_mapfile1_otutable1[complete.cases(reorder_mapfile1_otutable1),]");
    
    $R-> send("reorder_mapfile2_otutable2 <- mapfile2_otutable2[match(row.names(log_otutable2_L4),row.names(mapfile2_otutable2)),]");
    $R-> send("nona_reorder_mapfile2_otutable2 <- reorder_mapfile2_otutable2[complete.cases(reorder_mapfile2_otutable2),]");
    
    my $Rvar;
    
    # Trim the OTU table mapfiles to one another
    
    $Rcmdset = "nona_reorder2_mapfile1_otutable1 <- nona_reorder_mapfile1_otutable1[match(nona_reorder_mapfile2_otutable2\$Replicate,nona_reorder_mapfile1_otutable1\$Replicate),]
    nona2_reorder2_mapfile1_otutable1  <- nona_reorder2_mapfile1_otutable1[complete.cases(nona_reorder2_mapfile1_otutable1),]
    nona_reorder2_mapfile2_otutable2 <- nona_reorder_mapfile2_otutable2[match(nona2_reorder2_mapfile1_otutable1\$Replicate,nona_reorder_mapfile2_otutable2\$Replicate),]
    nona2_reorder2_mapfile2_otutable2 <- nona_reorder2_mapfile2_otutable2[complete.cases(nona_reorder2_mapfile2_otutable2),]";
    
    $R-> send("$Rcmdset");
        
    
    for(my $level=4;$level<8;$level++){
        
        
        # trim the corresponding OTU table to the OTUtable mapfile (which has just been re-trimmed).  Now both OTU tables rows should match.
        
        $Rcmdset = "reorder_log_otutable1_L".$level." <- log_otutable1_L".$level."[match(row.names(nona2_reorder2_mapfile1_otutable1), row.names(log_otutable1_L".$level.")),]
        nona_reorder_log_otutable1_L".$level." <- reorder_log_otutable1_L".$level."[complete.cases(reorder_log_otutable1_L".$level."),]
        reorder_otutable1_L".$level." <- otutable1_L".$level."[match(row.names(nona2_reorder2_mapfile1_otutable1), row.names(otutable1_L".$level.")),]
        nona_reorder_otutable1_L".$level." <- reorder_otutable1_L".$level."[complete.cases(reorder_otutable1_L".$level."),]";
        
        $R-> send("$Rcmdset");
        
        $Rcmdset = "reorder_log_otutable2_L".$level." <- log_otutable2_L".$level."[match(row.names(nona2_reorder2_mapfile2_otutable2), row.names(log_otutable2_L".$level.")),]
        nona_reorder_log_otutable2_L".$level." <- reorder_log_otutable2_L".$level."[complete.cases(reorder_log_otutable2_L".$level."),]
        reorder_otutable2_L".$level." <- otutable2_L".$level."[match(row.names(nona2_reorder2_mapfile2_otutable2), row.names(otutable2_L".$level.")),]
        nona_reorder_otutable2_L".$level." <- reorder_otutable2_L".$level."[complete.cases(reorder_otutable2_L".$level."),]";
        
        $R-> send("$Rcmdset");
    
    
        # Get the tübingen subset of each OTU table
        $Rcmdset = "tub_otutable1 <- nona_reorder_log_otutable1_L".$level."[which(nona2_reorder2_mapfile1_otutable1\$Tub_Loc != \"C\"),]
        tub_otutable2 <- nona_reorder_log_otutable2_L".$level."[which(nona2_reorder2_mapfile2_otutable2\$Tub_Loc != \"C\"),]";
        
        $R-> send("$Rcmdset");
        
        # Here with fixing up for dual OTU tables
        
        #$Rcmdset = "write.table(log_otutable1_L5, file=\"test.txt\", sep=\"\t\", quote = F)";
        #print $Rcmdset."\n";
        
        
        $nextRcmd = "rownames <- data.frame()
        rm(results_matrices_sums)
        rm(results_matrices_max)
        rm(results_matrices_numpos)
        rm(results_matrices_slopes)
        rm(results_matrices_rsq)
        rm(results_matrices_rsq_adj)
        rm(results_matrices_pvals)
        rm(results_matrices_pvals_adj)
        rm(results_matrices_tubsupport)
        rm(results_matrices_subsupport)
        
        for(i in 1:length(colnames(nona_reorder_log_otutable1_L".$level."))){
            
            otuvector1 <- nona_reorder_log_otutable1_L".$level."[,i]
            tub_otuvector1 <- tub_otutable1[,i]
            otusum1 <- sum(nona_reorder_otutable1_L".$level."[,i])
            otumax1 <- max(nona_reorder_otutable1_L".$level."[,i])
            smr1 <- otusum1 / otumax1
            numpos1 <-  length(which(otuvector1 > 0))
            
            if(otusum1 > 50 & numpos1 > 10){
                rownames <- c(rownames, colnames(nona_reorder_log_otutable1_L".$level.")[i])
                sums <- data.frame()
                maxes <- data.frame()
                numpos <- data.frame()
                slopes <- data.frame()
                rsq <- data.frame()
                rsq_adj <- data.frame()
                pvals <- data.frame()
                pvals_adj <- data.frame()
                tub_support <- data.frame()
                colnames <- data.frame()
                subsupport <- data.frame()
                
                for(j in 1:length(colnames(nona_reorder_log_otutable2_L".$level."))){
                    
                    otuvector2 <- nona_reorder_log_otutable2_L".$level."[,j]
                    tub_otuvector2 <- tub_otutable2[,j]
                    otusum2 <- sum(nona_reorder_otutable2_L".$level."[,j])
                    otumax2 <- max(nona_reorder_otutable2_L".$level."[,j])
                    smr2 <- otusum2 / otumax2
                    numpos2 <-  length(which(otuvector2 > 0))
                    
                    
                    
                    if(otusum2 > 50 & numpos2 > 10){
                        
                        resultlm <- summary(lm(otuvector1 ~ otuvector2))
                        slope <- resultlm\$coefficients[2,1]
                        r_square <- resultlm\$r.squared
                        adj_r_square <- resultlm\$adj.r.squared
                        if(slope < 0){
                            r_square <- r_square * -1
                            adj_r_square <- adj_r_square * -1
                        }
                        p_value <- pf(resultlm\$fstatistic[1], resultlm\$fstatistic[2], resultlm\$fstatistic[3], lower.tail=F)
                        p_value_adj <- p.adjust(p_value, method=\"fdr\",n=length(colnames(nona_reorder_log_otutable1_L".$level.")))
                        
                        
                        
                        resultlm_tub <- summary(lm(tub_otuvector1 ~ tub_otuvector2))
                        slope_tub <- resultlm_tub\$coefficients[2,1]
                        r_square_tub <- resultlm_tub\$r.squared
                        adj_r_square_tub <- resultlm_tub\$adj.r.squared
                        if(slope_tub < 0){
                            r_square_tub <- r_square_tub * -1
                            adj_r_square_tub <- adj_r_square_tub * -1
                        }
                        p_value_tub <- pf(resultlm_tub\$fstatistic[1], resultlm_tub\$fstatistic[2], resultlm_tub\$fstatistic[3], lower.tail=F)
                        p_value_adj_tub <- p.adjust(p_value_tub, method=\"fdr\",n=length(colnames(nona_reorder_log_otutable1_L".$level.")))
                        
                        if(p_value_tub < 0.05){
                            support = \"Y\"
                        }else{ support = \"N\" }
                        
                        write.table(c(i,j,p_value_adj), file=\"test0.txt\", sep=\"\t\", quote = F, append=TRUE)
                        
                        
                        sums <- c(sums, paste(otusum1,otusum2,sep=\"_\"))
                        maxes <- c(maxes, paste(otumax1,otumax2,sep=\"_\"))
                        numpos <- c(numpos,paste(numpos1,numpos2,sep=\"_\"))
                        slopes <- as.numeric(c(slopes, slope))
                        rsq <- as.numeric(c(rsq, r_square))
                        rsq_adj <- as.numeric(c(rsq_adj, adj_r_square))
                        pvals <- as.numeric(c(pvals, p_value))
                        pvals_adj <- as.numeric(c(pvals_adj, p_value_adj))
                        tub_support <- c(tub_support, support)
                        colnames <- c(colnames, colnames(nona_reorder_log_otutable2_L".$level.")[j])
                    }
                    
                }
                
                if(exists(\"results_matrices_sums\") != TRUE){
                    results_matrices_sums <- t(data.frame(1:length(sums)))
                    results_matrices_max <- t(data.frame(1:length(maxes)))
                    results_matrices_numpos <- t(data.frame(1:length(numpos)))
                    results_matrices_slopes <- t(data.frame(1:length(slopes)))
                    results_matrices_rsq <- t(data.frame(1:length(rsq)))
                    results_matrices_rsq_adj <- t(data.frame(1:length(rsq_adj)))
                    results_matrices_pvals <- t(data.frame(1:length(pvals)))
                    results_matrices_pvals_adj <- t(data.frame(1:length(pvals_adj)))
                    results_matrices_tubsupport <- t(data.frame(1:length(tub_support)))
                }
                
                results_matrices_sums <- rbind(results_matrices_sums, sums)
                results_matrices_max <- rbind(results_matrices_max, maxes)
                results_matrices_numpos <- rbind(results_matrices_numpos, numpos)
                results_matrices_slopes <- rbind(results_matrices_slopes, slopes)
                results_matrices_rsq <- rbind(results_matrices_rsq, rsq)
                results_matrices_rsq_adj <- rbind(results_matrices_rsq_adj, rsq_adj)
                results_matrices_pvals <- rbind(results_matrices_pvals, pvals)
                results_matrices_pvals_adj <- rbind(results_matrices_pvals_adj, pvals_adj)
                results_matrices_tubsupport <- rbind(results_matrices_tubsupport, tub_support)
            }
            
        }
        
        
        if(exists(\"results_matrices_sums\") == TRUE){
        
        row.names(results_matrices_sums) <- paste(row.names(results_matrices_sums), 1:length(row.names(results_matrices_sums)), sep=\"_\")
        row.names(results_matrices_max) <- paste(row.names(results_matrices_max), 1:length(row.names(results_matrices_max)), sep=\"_\")
        row.names(results_matrices_numpos) <- paste(row.names(results_matrices_numpos), 1:length(row.names(results_matrices_numpos)), sep=\"_\")
        row.names(results_matrices_slopes) <- paste(row.names(results_matrices_slopes), 1:length(row.names(results_matrices_slopes)), sep=\"_\")
        row.names(results_matrices_rsq) <- paste(row.names(results_matrices_rsq), 1:length(row.names(results_matrices_rsq)), sep=\"_\")
        row.names(results_matrices_rsq_adj) <- paste(row.names(results_matrices_rsq_adj), 1:length(row.names(results_matrices_rsq_adj)), sep=\"_\")
        row.names(results_matrices_pvals) <- paste(row.names(results_matrices_pvals), 1:length(row.names(results_matrices_pvals)), sep=\"_\")
        row.names(results_matrices_pvals_adj) <- paste(row.names(results_matrices_pvals_adj), 1:length(row.names(results_matrices_pvals_adj)), sep=\"_\")
        row.names(results_matrices_tubsupport) <- paste(row.names(results_matrices_tubsupport), 1:length(row.names(results_matrices_tubsupport)), sep=\"_\")
        
        if(nrow(results_matrices_sums) != 1){
            if(ncol(results_matrices_sums) == 1 & is.matrix(results_matrices_sums) == TRUE){
                results_matrices_sums <- data.frame(t(data.frame(results_matrices_sums[-1,])))
                results_matrices_max <- data.frame(t(data.frame(results_matrices_max[-1,])))
                results_matrices_numpos <- data.frame(t(data.frame(results_matrices_numpos[-1,])))
                results_matrices_slopes <- data.frame(results_matrices_slopes[-1,])
                results_matrices_rsq <- data.frame(results_matrices_rsq[-1,])
                results_matrices_rsq_adj <- data.frame(results_matrices_rsq_adj[-1,])
                results_matrices_pvals <- data.frame(results_matrices_pvals[-1,])
                results_matrices_pvals_adj <- data.frame(results_matrices_pvals_adj[-1,])
                results_matrices_tubsupport <- data.frame(t(data.frame(results_matrices_tubsupport[-1,])))
            
            }else if(nrow(results_matrices_sums) == 2){
                results_matrices_sums <- data.frame(results_matrices_sums[-1,])
                results_matrices_max <- data.frame(results_matrices_max[-1,])
                results_matrices_numpos <- data.frame(results_matrices_numpos[-1,])
                results_matrices_slopes <- data.frame(t(data.frame(results_matrices_slopes[-1,])))
                results_matrices_rsq <- data.frame(t(data.frame(results_matrices_rsq[-1,])))
                results_matrices_rsq_adj <- data.frame(t(data.frame(results_matrices_rsq_adj[-1,])))
                results_matrices_pvals <- data.frame(t(data.frame(results_matrices_pvals[-1,])))
                results_matrices_pvals_adj <- data.frame(t(data.frame(results_matrices_pvals_adj[-1,])))
                results_matrices_tubsupport <- data.frame(results_matrices_tubsupport[-1,])
            }
            else{
            results_matrices_sums <- results_matrices_sums[-1,]
            results_matrices_max <- results_matrices_max[-1,]
            results_matrices_numpos <- results_matrices_numpos[-1,]
            results_matrices_slopes <- results_matrices_slopes[-1,]
            results_matrices_rsq <- results_matrices_rsq[-1,]
            results_matrices_rsq_adj <- results_matrices_rsq_adj[-1,]
            results_matrices_pvals <- results_matrices_pvals[-1,]
            results_matrices_pvals_adj <- results_matrices_pvals_adj[-1,]
            results_matrices_tubsupport <- results_matrices_tubsupport[-1,]
            }
            
            results_matrices_sums[is.na(results_matrices_sums)] <- 0
            results_matrices_max[is.na(results_matrices_max)] <- 0
            results_matrices_numpos[is.na(results_matrices_numpos)] <- 0
            results_matrices_slopes[is.na(results_matrices_slopes)] <- 0
            results_matrices_rsq[is.na(results_matrices_rsq)] <- 0
            results_matrices_rsq_adj[is.na(results_matrices_rsq_adj)] <- 0
            results_matrices_pvals[is.na(results_matrices_pvals)] <- 0
            results_matrices_pvals_adj[is.na(results_matrices_pvals_adj)] <- 0
            results_matrices_tubsupport[is.na(results_matrices_tubsupport)] <- 0
            
            
            q <- c(colnames, ".$level.")
            colnames(results_matrices_sums) <- colnames
            colnames(results_matrices_max) <- colnames
            colnames(results_matrices_numpos) <- colnames
            colnames(results_matrices_slopes) <- colnames
            colnames(results_matrices_rsq) <- colnames
            colnames(results_matrices_rsq_adj) <- colnames
            colnames(results_matrices_pvals) <- colnames
            colnames(results_matrices_pvals_adj) <- colnames
            colnames(results_matrices_tubsupport) <- colnames
            
            
            
            row.names(results_matrices_sums) <- rownames
            row.names(results_matrices_max) <- rownames
            row.names(results_matrices_numpos) <- rownames
            row.names(results_matrices_slopes) <- rownames
            row.names(results_matrices_rsq) <- rownames
            row.names(results_matrices_rsq_adj) <- rownames
            row.names(results_matrices_pvals) <- rownames
            row.names(results_matrices_pvals_adj) <- rownames
            row.names(results_matrices_tubsupport) <- rownames
        }
        
        }";
        
        $R-> send("$nextRcmd"); 
        
        # Changed p-value filter to 1 so it doesn't filter anything anymore.
        
        $R-> send("if(exists(\"results_matrices_sums\") == TRUE){
        
        
            for(i in 1:length(colnames(results_matrices_sums))){
                for(j in 1:length(row.names(results_matrices_sums))){
                    
                    if(results_matrices_pvals[j,i] < 1){
                        
                        edge_vector <- data.frame()
                        
                        sumstr <- results_matrices_sums[j,i]
                        sum1 <- as.numeric(levels(as.data.frame(strsplit(as.character(sumstr),\"_\"))[,1]))[as.data.frame(strsplit(as.character(sumstr),\"_\"))[,1]][1]
                        sum2 <- as.numeric(levels(as.data.frame(strsplit(as.character(sumstr),\"_\"))[,1]))[as.data.frame(strsplit(as.character(sumstr),\"_\"))[,1]][2]
                        maxstr <- results_matrices_max[j,i]
                        max1 <- as.numeric(levels(as.data.frame(strsplit(as.character(maxstr),\"_\"))[,1]))[as.data.frame(strsplit(as.character(maxstr),\"_\"))[,1]][1]
                        max2 <- as.numeric(levels(as.data.frame(strsplit(as.character(maxstr),\"_\"))[,1]))[as.data.frame(strsplit(as.character(maxstr),\"_\"))[,1]][2]
                        numstr <- results_matrices_numpos[j,i]
                        num1 <- as.numeric(levels(as.data.frame(strsplit(as.character(numstr),\"_\"))[,1]))[as.data.frame(strsplit(as.character(numstr),\"_\"))[,1]][1]
                        num2 <- as.numeric(levels(as.data.frame(strsplit(as.character(numstr),\"_\"))[,1]))[as.data.frame(strsplit(as.character(numstr),\"_\"))[,1]][2]
                        
                        stm1 <- sum1 / max1
                        stm2 <- sum2 / max2
                        avestm <- (stm1 + stm2) / 2
                        stm_rsqadj <- abs(avestm * results_matrices_rsq_adj[j,i])
                        
                        org1 <- row.names(results_matrices_sums)[j]
                        org1 <- paste(\"".$endoepi1."\",org1,sep = \"_\")
                        org2 <- colnames(results_matrices_sums)[i]
                        org2 <- paste(\"".$endoepi2."\",org2,sep = \"_\")
                        
                        edge_vector <- c(org1,org2,results_matrices_tubsupport[j,i],\"".$groupname1."_".$groupname2."\",results_matrices_slopes[j,i],results_matrices_rsq[j,i],results_matrices_rsq_adj[j,i],results_matrices_pvals[j,i],results_matrices_pvals_adj[j,i],avestm,stm_rsqadj)
                        edgetable_".$level." <- rbind(edgetable_".$level.", edge_vector)
                        
                        if(TRUE %in% grepl(org1, nodetable_".$level."[,1])){
                        }else{
                        orgsplit <- as.data.frame(strsplit(org1, \";\"))
                        kingdom <- as.character(orgsplit[1,])
                        phyla <- as.character(orgsplit[2,])
                        class <- as.character(orgsplit[3,])
                        order <- as.character(orgsplit[4,])
                        if(".$level." > 4){family <- as.character(orgsplit[5,])}else{family <- \"n/a\"}
                        if(".$level." > 5){genus <- as.character(orgsplit[6,])}else{genus <- \"n/a\"}
                        orgnodeaddition1 <- data.frame()
                        orgnodeaddition1 <- c(org1,\"Organism\",\"".$endoepi1."\",\"".$groupname1."\", kingdom, phyla, class, order, family, genus, sum1, max1, num1)
                        nodetable_".$level." <- rbind(nodetable_".$level.", orgnodeaddition1)
                        }
                        
                        if(TRUE %in% grepl(org2, nodetable_".$level."[,1])){
                        }else{
                        orgsplit <- as.data.frame(strsplit(org2, \";\"))
                        kingdom <- as.character(orgsplit[1,])
                        phyla <- as.character(orgsplit[2,])
                        class <- as.character(orgsplit[3,])
                        order <- as.character(orgsplit[4,])
                        if(".$level." > 4){family <- as.character(orgsplit[5,])}else{family <- \"n/a\"}
                        if(".$level." > 5){genus <- as.character(orgsplit[6,])}else{genus <- \"n/a\"}
                        orgnodeaddition1 <- data.frame()
                        orgnodeaddition1 <- c(org2,\"Organism\",\"".$endoepi2."\",\"".$groupname2."\", kingdom, phyla, class, order, family, genus, sum2, max2, num2)
                        nodetable_".$level." <- rbind(nodetable_".$level.", orgnodeaddition1)
                        }
                        
                    }
                    
                }
                
            }
        
        }");
        
    }
    
    for(my $i=4;$i<8;$i++){
    
    my $Rcmd = "if(exists(\"results_matrices_sums\") == TRUE){
    
    write.table(edgetable_".$i.", file=\"".$resultsdir."EdgeTable".$i.".txt\", sep=\"\t\", quote=F, row.names=F)
    write.table(nodetable_".$i.", file=\"".$resultsdir."NodeTable".$i.".txt\", sep=\"\t\", quote=F, row.names=F)
    
    }";
    
    $R-> send("$Rcmd");
    
    }   
    
    $R->stopR();
    
}

exit;

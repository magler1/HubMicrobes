# First run export PERL5LIB=$PERL5LIB:/Users/agler/packages/Statistics-R-0.33/lib/:/Users/agler/packages/IPC-Run-0.94/lib/:/Users/agler/packages/Regexp-Common-2013031301/lib/
# export R_LIBS="/Users/agler/R-packages/"

# Resubmission Version: This version filters out the May 5 time point from the cologne samples.

use Statistics::R;

my $otutablelist = shift; # File listing tables to summarize and do ordination on, and tab deliminating the depth at which to filter that table

open TABLELIST, $otutablelist;

while(<TABLELIST>){

    my $currentline = $_;
    chomp $currentline;
    my @splitline = split(/\t/,$currentline);
    my $otutable = "OTU_Tables/".$splitline[0];
    my $depth = $splitline[1];
    my $mapfile = "Mapfiles/".$splitline[2];
    my @splitmapfile = split(/\./,$mapfile);
    my $mapfilebasename = $splitmapfile[0];
    
    my @tablenamesplit = split(/\./,$otutable);
    my $tablebasename = $tablenamesplit[0];
    my @basenamesplit = split(/\//,$tablebasename);
    my $realbase = $basenamesplit[1];
    
    my $rarefiedtable = $tablebasename."_".$depth.".biom";
    my $lookfor1 = $rarefiedtable;
    
    my $command1 = "single_rarefaction.py -i $otutable -o $rarefiedtable -d $depth";
    my $cmd1 = system('bash','-c',". /macqiime/configs/bash_profile.txt && $command1") == 0
    or die "system failed: $?";
    sleep 1 until -e $lookfor1;
    
    my $summarydir = $tablebasename."_".$depth."_taxsumm/";
    my $lookfor2 = $summarydir.$realbase."_".$depth."_L6.txt";
    
    my $command2 = "summarize_taxa.py -i $rarefiedtable -o $summarydir -a";
    my $cmd2 = system('bash', '-c', ". /macqiime/configs/bash_profile.txt && $command2") == 0
    or die "system failed: $?";
    sleep 1 until -e $lookfor2;

    
    # Trim mapfile and put it samples in the same order as the OTU table
    
    #First load table into multi-dim hash of form $hash{$samplename}{$otuname}=$obscount;
    
    my $summtable = $summarydir.$realbase."_".$depth."_L4.txt";         # It doesn't matter which table to use, samples are all trimmed in the same way

    open TEXTTABLE, $summtable;
    
    my @sampleorder;
    
    while(<TEXTTABLE>){
        
        my $currentline = $_;
        chomp $currentline;
        my @splitline = split(/\t/,$currentline);
        
        
        if($currentline =~ /Taxon/){
            
            $j=0;
            foreach(@splitline){
                if($j==0){$j++; next;}
                my $q = $j - 1;
             $sampleorder[$q] = $splitline[$j];
             $j++   
            }
        }
        
        next;
        
    }
    close(TEXTTABLE);
    
    ### Now open the mapfile and re-order the samples, get rid of lines that are not in otu table
    
    my $trimmedmapfile = $mapfilebasename."_trimmap.txt";
    open MAPFILE, $mapfile;
    open MAPFILE2, ">".$trimmedmapfile;
    
    my %mapfilehash;
    
    while(<MAPFILE>){
        
        my $currentline = $_;
        chomp $currentline;

        if ($currentline =~ /#/){
            print MAPFILE2 "$currentline\n";
            next;
        }
        
        # Load the mapfile into a hash with sample names as keys
        my @splitline = split(/\t/,$currentline);
    
        $mapfilehash{$splitline[0]} = $currentline;
        
    }
    
    close(MAPFILE);
    
    foreach(@sampleorder){

        my $printline = $mapfilehash{$_};
        print MAPFILE2 "$printline\n";
        
    }
    
    close(MAPFILE2);
    
    my $resultsdir = $summarydir."Results/";
    if(-e $resultsdir){}else{
        my $command3 = "mkdir $resultsdir";
        my $cmd3 = system('bash', '-c', "$command3") == 0
        or die "system failed: $?";
    }
    my $colresultsdir = $resultsdir."Cologne_only/";
    if(-e $colresultsdir){}else{
        my $command4 = "mkdir $colresultsdir";
        my $cmd4 = system('bash', '-c', "$command4") == 0
        or die "system failed: $?";
    }
    
    for(my $i=4;$i<7;$i++){
        
        my $R = Statistics::R->new();
    
        $R->startR;
        
        my $currenttable = $summarydir.$realbase."_".$depth."_L".$i.".txt";
        
        my $outresults = $resultsdir.$realbase."_".$depth."_L".$i."_CA_sigResults.txt";
        my $outcoords = $resultsdir.$realbase."_".$depth."_L".$i."_CA_coordseigs.txt";
        my $logoutresults = $resultsdir.$realbase."_".$depth."_L".$i."_CA_sigResults_log.txt";
        my $logoutcoords = $resultsdir.$realbase."_".$depth."_L".$i."_CA_coordseigs_log.txt";
        
        $R-> send("otutable <- read.table(\"$currenttable\", row.names=1, comment.char=\"\", sep=\"\t\", header=T)");
        $R-> send(q`otutable <- t(otutable)`);
        
        
        $R-> send("mapfile <- read.table(\"$trimmedmapfile\",header=T, row.names=1, comment.char=\"\", sep=\"\t\")");   # Feed in the trimmed mapfiles and perform CCA to get strength of effects of factors
    
        # Remove the May 5 data points in the cologne data set, since we don't have even sampling from this timepoint
    
        $R-> send(q`otutable <- subset(otutable, mapfile$Tub_Loc_Time!="5")`);
        $R-> send(q`log_otutable <- log10(otutable + 1)`);
        $R-> send(q`mapfile <- subset(mapfile, mapfile$Tub_Loc_Time!="5")`);
    
        $R-> send(q`library(vegan)`);
        $R-> send(q`library(gdata)`);
        $R-> send(q`library(plyr)`);
        $R-> send(q`library(ggplot2)`);
        $R-> send(q`library(gridExtra)`);
        $R-> send(q`library(colorspace)`);
        $R-> send(q`library(RColorBrewer)`);
       
        
        # Also cologne only by ecotype and by samptime would be nice...
        
        $R-> send(q`otutable_col <- subset(otutable, mapfile$Lab_Col_Tub=="C")`);
        $R-> send(q`log_otutable_col <- log10(otutable_col + 1)`);
        $R-> send(q`mapfile_col <- subset(mapfile, mapfile$Lab_Col_Tub=="C")`);
        
         
        
            # Cologne only unconstrained
        $R-> send("otutable_col_ca <- cca(otutable_col)");                           # Run the cca on cologne samples only
        $R-> send("log_otutable_col_ca <- cca(log_otutable_col)");
        $R-> send(q`CA1perc_col <- (otutable_col_ca$CA$eig[1]/sum(otutable_col_ca$CA$eig))*100`);
        $R-> send(q`CA1perc_col <- round(CA1perc_col, digits = 1)`);
        $R-> send(q`CA2perc_col <- (otutable_col_ca$CA$eig[2]/sum(otutable_col_ca$CA$eig))*100`);
        $R-> send(q`CA2perc_col <- round(CA2perc_col, digits = 1)`);
        $R-> send(q`CA1perc_col_log <- (log_otutable_col_ca$CA$eig[1]/sum(log_otutable_col_ca$CA$eig))*100`);
        $R-> send(q`CA1perc_col_log <- round(CA1perc_col_log, digits = 1)`);
        $R-> send(q`CA2perc_col_log <- (log_otutable_col_ca$CA$eig[2]/sum(log_otutable_col_ca$CA$eig))*100`);
        $R-> send(q`CA2perc_col_log <- round(CA2perc_col_log, digits = 1)`);
        
            # Cologne only Ecotype Constrained
        $R-> send("otutable_col_cca_Ecotype <- cca(otutable_col ~ Ecotype, data=mapfile_col)");
       
        $R-> send(q`CCA1perc_Ecotype <- (sum(otutable_col_cca_Ecotype$CCA$eig[1])/sum(otutable_col_cca_Ecotype$CCA$eig,otutable_col_cca_Ecotype$CA$eig))*100`);
        $R-> send(q`CCA1perc_Ecotype <- round(CCA1perc_Ecotype, digits = 1)`);
        $R-> send(q`CCA2perc_Ecotype <- (sum(otutable_col_cca_Ecotype$CCA$eig[2])/sum(otutable_col_cca_Ecotype$CCA$eig,otutable_col_cca_Ecotype$CA$eig))*100`);
        $R-> send(q`CCA2perc_Ecotype <- round(CCA2perc_Ecotype, digits = 1)`);
        $R-> send(q`CCAperc_Ecotype <- (sum(otutable_col_cca_Ecotype$CCA$eig)/sum(otutable_col_cca_Ecotype$CCA$eig,otutable_col_cca_Ecotype$CA$eig))*100`);  # Run the cca
        $R-> send(q`CCAperc_Ecotype <- round(CCAperc_Ecotype, digits = 1)`);
        $R-> send("log_otutable_col_cca_Ecotype <- cca(log_otutable_col ~ Ecotype, data=mapfile_col)");
        $R-> send(q`CCAperc_log_Ecotype <- (sum(log_otutable_col_cca_Ecotype$CCA$eig)/sum(log_otutable_col_cca_Ecotype$CCA$eig,log_otutable_col_cca_Ecotype$CA$eig))*100`);
        $R-> send(q`CCAperc_log_Ecotype <- round(CCAperc_log_Ecotype, digits = 1)`);
        $R-> send(q`CA1perc_col_log_Ecotype <- (log_otutable_col_cca_Ecotype$CA$eig[1]/sum(log_otutable_col_cca_Ecotype$CA$eig))*100`);
        $R-> send(q`CA1perc_col_log_Ecotype <- round(CA1perc_col_log_Ecotype, digits = 1)`);
        $R-> send(q`CCA1perc_log_Ecotype <- (sum(log_otutable_col_cca_Ecotype$CCA$eig[1])/sum(log_otutable_col_cca_Ecotype$CCA$eig,log_otutable_col_cca_Ecotype$CA$eig))*100`);
        $R-> send(q`CCA1perc_log_Ecotype <- round(CCA1perc_log_Ecotype, digits = 1)`);
        $R-> send(q`CCA2perc_log_Ecotype <- (sum(log_otutable_col_cca_Ecotype$CCA$eig[2])/sum(log_otutable_col_cca_Ecotype$CCA$eig,log_otutable_col_cca_Ecotype$CA$eig))*100`);
        $R-> send(q`CCA2perc_log_Ecotype <- round(CCA2perc_log_Ecotype, digits = 1)`);
        
        $nextoutfile = $colresultsdir.$realbase."_".$depth."_L".$i."_Ecotype_SigTest.txt";
        
        $R-> send("SigTest1 <- anova(log_otutable_col_cca_Ecotype, by=\"term\", step=9999, perm.max=9999)");
        $R-> send("SigTest1_axes <- anova(log_otutable_col_cca_Ecotype, by=\"axis\", step=9999, perm.max=9999)");
        
        $R-> send("LocSTResults <- c(CCAperc_log_Ecotype, SigTest1\$\'Pr(>F)\'[1], SigTest1_axes\$\'Pr(>F)\'[1], SigTest1_axes\$\'Pr(>F)\'[2])");
        
        $R-> send("Results <- rbind(LocSTResults)
                  colnames(Results) <- c(\"Percent_Correlated\", \"pval_Ecotype\", \"pval_CCA1\", \"pval_CCA2\")");
        
        $R-> send("write.table(Results, file=\"".$nextoutfile."\", sep=\"\t\", quote=F, col.names=NA)");
        
        # Set up the themes for the plots
        $R-> send("t2 <- theme(
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  
                  legend.key = element_blank(),
                  legend.position = \"right\",
                  
                  axis.ticks.length = unit(-0.25, \"cm\"),
                  axis.ticks.margin = unit(0.5, \"cm\"),
                  axis.ticks = element_line(colour= \"black\", size=2),
                  axis.text = element_blank(),
                  axis.line = element_line(colour = \"black\", size = 2),
                  axis.title= element_text(face=\"bold\", color=\"black\", size=17),
                  
                  plot.title = element_text(face=\"bold\", color = \"black\", size=17))");
        
        $R-> send("getPalette = colorRampPalette(brewer.pal(8, \"Accent\"))
                    myColors <- getPalette(6)
                    myColors2 <- getPalette(3)
                    Loc <- c(\"PFN\",\"ERG\",\"WH\",\"EY\",\"JUG\",\"C\")
                    Eco <- c(\"C\",\"K\",\"W\")
                    df <- data.frame(1:6)
                    df2 <- data.frame(1:3)
                    df\$Loc <- as.factor(Loc)
                    df2\$Eco <- as.factor(Eco)
                    names(myColors) <- levels(df\$Loc)
                    names(myColors2) <- levels(df2\$Eco)");
        
        
        # Make the Cologne unconstrained plots
        
        $R-> send(q`data_for_plot <- cbind(log_otutable_col_ca$CA$u,mapfile_col)`);
        
        my $outfilename1 = $colresultsdir.$realbase."_dataforplot_L".$i.".txt";
        $R-> send("write.table(data_for_plot, file=\"$outfilename1\", sep=\"\t\", quote=F, col.names=NA)");
        
        $R-> send(q`gp_uncon_log_Col_Eco_jitt <- ggplot() + geom_point(data=data_for_plot, aes(x=CA1, y=CA2, color=Ecotype, shape=SampTime), size=7, position = position_jitter(w = 0.1, h = 0.1))`);
        $R-> send("gp_uncon_log_Col_Eco_jitt_adj <- gp_uncon_log_Col_Eco_jitt + t2 + scale_colour_manual(name=\"Eco\", values=myColors2) + labs(x=paste(\"CA1 (\",CA1perc_col_log,\"% Variation)\"), y =paste(\"CA2 (\",CA2perc_col_log,\"% Variation)\"))");
        
        $R-> send(q`gp_uncon_log_Col_Eco <- ggplot() + geom_point(data=data_for_plot, aes(x=CA1, y=CA2, color=Ecotype, shape=SampTime), size=7, position = position_jitter(w = 0.1, h = 0.1))`);
        $R-> send("gp_uncon_log_Col_Eco_adj <- gp_uncon_log_Col_Eco + t2 + scale_colour_manual(name=\"Eco\", values=myColors2) + labs(x=paste(\"CA1 (\",CA1perc_col_log,\"% Variation)\"), y =paste(\"CA2 (\",CA2perc_col_log,\"% Variation)\"))");
        
        $R-> send(q`gp_uncon_log_Col_Strain <- ggplot() + geom_point(data=data_for_plot, aes(x=CA1, y=CA2, color=Albugo_Strain, shape=SampTime), size=7, position = position_jitter(w = 0.1, h = 0.1))`);
        $R-> send("gp_uncon_log_Col_Strain_adj <- gp_uncon_log_Col_Strain + t2 + labs(x=paste(\"CA1 (\",CA1perc_col_log,\"% Variation)\"), y =paste(\"CA2 (\",CA2perc_col_log,\"% Variation)\"))");
        
        my $outfilename = $colresultsdir.$realbase."_Uncon_Col_log_L".$i.".pdf";
        $R-> send("pdf(\"$outfilename\", width=10, height=3.5)");
        $R-> send(q`grid.arrange(gp_uncon_log_Col_Eco_adj,gp_uncon_log_Col_Strain_adj,nrow=1)`);
        $R-> send("dev.off()");
        
        # Make the Cologne Ecotype constrained plots

        $R-> send(q`data_for_plot <- cbind(log_otutable_col_cca_Ecotype$CCA$wa,log_otutable_col_cca_Ecotype$CA$u,mapfile_col)`);
        
         
        $R-> send("if(length(colnames(log_otutable_col_cca_Ecotype\$CCA\$wa)) > 1){
                        gp_EcoCon_log_Col_Eco <- ggplot() + geom_point(data=data_for_plot, aes(x=CCA1, y=CCA2, color=Ecotype, shape=SampTime), size=7, position = position_jitter(w = 0.1, h = 0.1))
                        gp_EcoCon_log_Col_Eco_adj <- gp_EcoCon_log_Col_Eco + t2 + scale_colour_manual(name=\"Eco\", values=myColors2) + labs(title=paste(\"Total constrained: \",CCAperc_log_Ecotype,\"% Variation\"),x=paste(\"CCA1 (\",CCA1perc_log_Ecotype,\"% Variation)\"), y =paste(\"CCA2 (\",CCA2perc_log_Ecotype,\"% Variation)\")) + theme(plot.title = element_text(vjust=3))
                    }else{
                        gp_EcoCon_log_Col_Eco <- ggplot() + geom_point(data=data_for_plot, aes(x=CCA1, y=CA1, color=Ecotype, shape=SampTime), size=7, position = position_jitter(w = 0.1, h = 0.1))
                        gp_EcoCon_log_Col_Eco_adj <- gp_EcoCon_log_Col_Eco + t2 + scale_colour_manual(name=\"Eco\", values=myColors2) + labs(title=paste(\"Total constrained: \",CCAperc_log_Ecotype,\"% Variation\"),x=paste(\"CCA1 (\",CCA1perc_log_Ecotype,\"% Variation)\"), y =paste(\"CA1 (\",CA1perc_col_log_Ecotype,\"% Variation)\")) + theme(plot.title = element_text(vjust=3))
                    }");
        
        my $outfilename = $colresultsdir.$realbase."_EcoCon_log_Col_L".$i.".pdf";
        $R-> send("pdf(\"$outfilename\", width=5, height=3.5)");
        $R-> send(q`grid.arrange(gp_EcoCon_log_Col_Eco_adj,nrow=1)`);
        $R-> send("dev.off()");
        
        $R->stopR();
    }
    
}

close TABLELIST;

exit;

# First run export PERL5LIB=$PERL5LIB:/Users/agler/packages/Statistics-R-0.33/lib/:/Users/agler/packages/IPC-Run-0.94/lib/:/Users/agler/packages/Regexp-Common-2013031301/lib/
# export R_LIBS="/Users/agler/R-packages/"

# Resubmission Version: This version filters out the May 5 time point from the cologne samples because these were not evenly sampled.

use Statistics::R;

my $otutablelist = shift; # File listing tables to summarize and do ordination on, and tab deliminating the depth at which to filter that table
my $plotlevel = shift; # 4,5, or 6 : level at which to make a plot of the principal coordinate plots
my $plotcat = shift; # Metadata category to use when coloring the plots

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
    
    if(-e $lookfor1){}else{
    my $command1 = "single_rarefaction.py -i $otutable -o $rarefiedtable -d $depth";
    my $cmd1 = system('bash','-c',". /macqiime/configs/bash_profile.txt && $command1") == 0
    or die "system failed: $?";
    sleep 1 until -e $lookfor1;
    }
    
    my $summarydir = $tablebasename."_".$depth."_taxsumm/";
    my $lookfor2 = $summarydir.$realbase."_".$depth."_L6.txt";
    
    if(-e $lookfor2){}else{
    my $command2 = "summarize_taxa.py -i $rarefiedtable -o $summarydir -a";
    my $cmd2 = system('bash', '-c', ". /macqiime/configs/bash_profile.txt && $command2") == 0
    or die "system failed: $?";
    sleep 1 until -e $lookfor2;
    }
    
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
    #print join(", ", @sampleorder);
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
        $R-> send(q`log_otutable <- log10(otutable + 1)`);
        
        $R-> send("mapfile <- read.table(\"$trimmedmapfile\",header=T, row.names=1, comment.char=\"\", sep=\"\t\")");   # Feed in the trimmed mapfiles and perform CCA to get strength of effects of factors
    
        # Since this is for tubingen samples only, subset the OTU table and mapfile to those datasets before doing the analysis
        $R-> send(q`otutable <- subset(otutable, mapfile$Lab_Col_Tub=="T")`);
        $R-> send(q`log_otutable <- log10(otutable + 1)`);
        $R-> send(q`mapfile <- subset(mapfile, mapfile$Lab_Col_Tub=="T")`);
    
        $R-> send(q`library(vegan)`);
        $R-> send(q`library(plyr)`);
        $R-> send(q`library(ggplot2)`);
        $R-> send(q`library(gridExtra)`);
        $R-> send(q`library(colorspace)`);
        $R-> send(q`library(RColorBrewer)`);
    
        $R-> send("otutable_ca <- cca(otutable)");                           # Run the cca
        $R-> send("log_otutable_ca <- cca(log_otutable)");
        
        $R-> send(q`coords <- otutable_ca$CA$u`);                                 # get the coordinates of samples and eigenvalues of vectors
        $R-> send(q`log_coords <- log_otutable_ca$CA$u`);
        $R-> send(q`eigvals <- otutable_ca$CA$eig`);
        $R-> send(q`CA1perc <- (otutable_ca$CA$eig[1]/sum(otutable_ca$CA$eig))*100`);
        $R-> send(q`CA1perc <- round(CA1perc, digits = 1)`);
        $R-> send(q`CA2perc <- (otutable_ca$CA$eig[2]/sum(otutable_ca$CA$eig))*100`);
        $R-> send(q`CA2perc <- round(CA2perc, digits = 1)`);
        $R-> send(q`log_eigvals <- log_otutable_ca$CA$eig`);
        $R-> send(q`CCA_data <- rbind(coords,eigvals)`);
        $R-> send(q`log_CCA_data <- rbind(log_coords,log_eigvals)`);
        $R-> send(q`CA1perc_log <- (log_otutable_ca$CA$eig[1]/sum(log_otutable_ca$CA$eig))*100`);
        $R-> send(q`CA1perc_log <- round(CA1perc_log, digits = 1)`);
        $R-> send(q`CA2perc_log <- (log_otutable_ca$CA$eig[2]/sum(log_otutable_ca$CA$eig))*100`);
        $R-> send(q`CA2perc_log <- round(CA2perc_log, digits = 1)`);
    
        $R-> send("write.table(log_CCA_data, file=\"$logoutcoords\", sep=\"\t\", quote=F, col.names=NA)");
        
        # Constrained by Tub_Loc (essentially location)
        $R-> send("otutable_cca_Loc <- cca(otutable ~ Tub_Loc, data=mapfile)");                           # Run the cca
        $R-> send(q`CCA1perc_Loc <- (sum(otutable_cca_Loc$CCA$eig[1])/sum(otutable_cca_Loc$CCA$eig,otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCA1perc_Loc <- round(CCA1perc_Loc, digits = 1)`);
        $R-> send(q`CCA2perc_Loc <- (sum(otutable_cca_Loc$CCA$eig[2])/sum(otutable_cca_Loc$CCA$eig,otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCA2perc_Loc <- round(CCA2perc_Loc, digits = 1)`);
        $R-> send(q`CCAperc_Loc <- (sum(otutable_cca_Loc$CCA$eig)/sum(otutable_cca_Loc$CCA$eig,otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCAperc_Loc <- round(CCAperc_Loc, digits = 1)`);
        $R-> send("log_otutable_cca_Loc <- cca(log_otutable ~ Tub_Loc, data=mapfile)");
        $R-> send(q`CCA1perc_log_Loc <- (sum(log_otutable_cca_Loc$CCA$eig[1])/sum(log_otutable_cca_Loc$CCA$eig,log_otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCA1perc_log_Loc <- round(CCA1perc_log_Loc, digits = 1)`);
        $R-> send(q`CCA2perc_log_Loc <- (sum(log_otutable_cca_Loc$CCA$eig[2])/sum(log_otutable_cca_Loc$CCA$eig,log_otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCA2perc_log_Loc <- round(CCA2perc_log_Loc, digits = 1)`);
        $R-> send(q`CCAperc_log_Loc <- (sum(log_otutable_cca_Loc$CCA$eig)/sum(log_otutable_cca_Loc$CCA$eig,log_otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCAperc_log_Loc <- round(CCAperc_log_Loc, digits = 1)`);
        
        $R-> send("t1 <- theme(
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  
                  legend.key = element_blank(),
                  legend.position = \"right\",
                  
                  axis.ticks.length = unit(-0.25, \"cm\"),
                  axis.ticks.margin = unit(0.5, \"cm\"),
                  axis.ticks.y = element_line(colour= \"black\", size=2),
                  axis.ticks.x = element_blank(),
                  axis.text = element_text(face=\"bold\", color=\"black\", size=25),
                  axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
                  axis.line = element_line(colour = \"black\", size = 2),
                  axis.title= element_text(face=\"bold\", color=\"black\", size=10),
                  
                  plot.title = element_blank())");
        
        # Constrained by SampTime (Fall/Spring and two times in Cologne)
        $R-> send("otutable_cca_SampTime <- cca(otutable ~ SampTime, data=mapfile)");
        $R-> send(q`CCA1perc_SampTime <- (sum(otutable_cca_SampTime$CCA$eig[1])/sum(otutable_cca_SampTime$CCA$eig,otutable_cca_SampTime$CA$eig))*100`);
        $R-> send(q`CCA1perc_SampTime <- round(CCA1perc_SampTime, digits = 1)`);
        $R-> send(q`CA1perc_SampTime <- (sum(otutable_cca_SampTime$CA$eig[1])/sum(otutable_cca_SampTime$CCA$eig,otutable_cca_SampTime$CA$eig))*100`);
        $R-> send(q`CA1perc_SampTime <- round(CA1perc_SampTime, digits = 1)`);
        $R-> send(q`CCAperc_SampTime <- (sum(otutable_cca_SampTime$CCA$eig)/sum(otutable_cca_SampTime$CCA$eig,otutable_cca_SampTime$CA$eig))*100`);  # Run the cca
        $R-> send(q`CCAperc_SampTime <- round(CCAperc_SampTime, digits = 1)`);
        $R-> send("log_otutable_cca_SampTime <- cca(log_otutable ~ SampTime, data=mapfile)");
        $R-> send(q`CCAperc_log_SampTime <- (sum(log_otutable_cca_SampTime$CCA$eig)/sum(log_otutable_cca_SampTime$CCA$eig,log_otutable_cca_SampTime$CA$eig))*100`);
        $R-> send(q`CCAperc_log_SampTime <- round(CCAperc_log_SampTime, digits = 1)`);
        $R-> send(q`CCA1perc_log_SampTime <- (sum(log_otutable_cca_SampTime$CCA$eig[1])/sum(log_otutable_cca_SampTime$CCA$eig,log_otutable_cca_SampTime$CA$eig))*100`);
        $R-> send(q`CCA1perc_log_SampTime <- round(CCA1perc_log_SampTime, digits = 1)`);
        $R-> send(q`CA1perc_log_SampTime <- (sum(log_otutable_cca_SampTime$CA$eig[1])/sum(log_otutable_cca_SampTime$CCA$eig,log_otutable_cca_SampTime$CA$eig))*100`);
        $R-> send(q`CA1perc_log_SampTime <- round(CA1perc_log_SampTime, digits = 1)`);
        
        # Constrained by Tub_Loc and SampTime (to check for independence)
        $R-> send("otutable_cca_Loc_ST <- cca(otutable ~ Tub_Loc + SampTime, data=mapfile)");
        $R-> send(q`CCAperc_Loc <- (sum(otutable_cca_Loc$CCA$eig)/sum(otutable_cca_Loc$CCA$eig,otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCAperc_Loc <- round(CCAperc_Loc, digits = 1)`);
        $R-> send("log_otutable_cca_Loc_ST <- cca(log_otutable ~ Tub_Loc + SampTime, data=mapfile)");
        $R-> send(q`CCAperc_log_Loc <- (sum(log_otutable_cca_Loc$CCA$eig)/sum(log_otutable_cca_Loc$CCA$eig,log_otutable_cca_Loc$CA$eig))*100`);
        $R-> send(q`CCAperc_log_Loc <- round(CCAperc_log_Loc, digits = 1)`);
        
        $nextoutfile = $resultsdir.$realbase."_".$depth."_L".$i."_LocationST_SigTest.txt";
        
        $R-> send("SigTest1 <- anova(log_otutable_cca_Loc_ST, by=\"term\", step=9999, perm.max=9999)");
        $R-> send("SigTest1_axis <- anova(log_otutable_cca_Loc_ST, by=\"axis\", step=9999, perm.max=9999)");
        $R-> send("SigTest2 <- anova(log_otutable_cca_Loc, by=\"term\", step=9999, perm.max=9999)");
        $R-> send("SigTest3 <- anova(log_otutable_cca_SampTime, by=\"term\", step=9999, perm.max=9999)");
        $R-> send("SigTest2_axis <- anova(log_otutable_cca_Loc, by=\"axis\", step=9999, perm.max=9999)");
        
        $R-> send(q`CCAperc_log_LocST <- (sum(log_otutable_cca_Loc_ST$CCA$eig)/sum(log_otutable_cca_Loc_ST$CCA$eig,log_otutable_cca_Loc_ST$CA$eig))*100`);
        $R-> send("LocSTResults <- c(CCAperc_log_LocST, SigTest1\$\'Pr(>F)\'[1], SigTest1\$\'Pr(>F)\'[2], SigTest1_axis\$\'Pr(>F)\'[1], SigTest1_axis\$\'Pr(>F)\'[2])");
        
        $R-> send(q`CCAperc_log_Loc <- (sum(log_otutable_cca_Loc$CCA$eig)/sum(log_otutable_cca_Loc$CCA$eig,log_otutable_cca_Loc$CA$eig))*100`);
        $R-> send("LocResults <- c(CCAperc_log_Loc, SigTest2\$\'Pr(>F)\'[1], \"NA\", SigTest2_axis\$\'Pr(>F)\'[1], SigTest2_axis\$\'Pr(>F)\'[2])");
        
        $R-> send(q`CCAperc_log_SampTime <- (sum(log_otutable_cca_SampTime$CCA$eig)/sum(log_otutable_cca_SampTime$CCA$eig,log_otutable_cca_SampTime$CA$eig))*100`);
        $R-> send("STResults <- c(CCAperc_log_SampTime, \"NA\", SigTest3\$\'Pr(>F)\'[1], \"NA\", \"NA\")");
        
        $R-> send("Results <- rbind(LocSTResults,LocResults,STResults)
                  colnames(Results) <- c(\"Percent_Correlated\", \"pval_Loc\", \"pval_SampTime\", \"CCA1_pval\", \"CCA2_pval\")");
        
        $R-> send("write.table(Results, file=\"".$nextoutfile."\", sep=\"\t\", quote=F, col.names=NA)");
        
        # PLot the principal coordinates for the current dataset for unconstrained_log (location), unconstrained_log (Albugo), and the constraints colored by location or the constraints themselves
        
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
        
            $R-> send("getPalette = colorRampPalette(brewer.pal(6, \"Dark2\"))
            myColors <- getPalette(6)
            Loc <- c(\"PFN\",\"ERG\",\"WH\",\"EY\",\"JUG\",\"C\")
            df <- data.frame(1:6)
            df\$Loc <- as.factor(Loc)
            names(myColors) <- levels(df\$Loc)");
        
        # Make the unconstrained plots
        $R-> send(q`data_for_plot <- cbind(log_otutable_ca$CA$u,mapfile)`);
        
        $R-> send(q`gp_uncon_log_Loc <- ggplot() + geom_point(data=data_for_plot, aes(x=CA1, y=CA2, color=Tub_Loc, shape=SampTime), size=7)`);
        $R-> send("gp_uncon_log_Loc_adj <- gp_uncon_log_Loc + t2 + scale_colour_manual(name=\"Loc\", values=myColors) + labs(x=paste(\"CA1 (\",CA1perc_log,\"% Variation)\"), y =paste(\"CA2 (\",CA2perc_log,\"% Variation)\"))");
        
        my $outfilename = $resultsdir.$realbase."_Uncon_log_L".$i.".pdf";
        $R-> send("pdf(\"$outfilename\", width=5, height=3.5)");
        $R-> send(q`grid.arrange(gp_uncon_log_Loc_adj,nrow=1)`);
        $R-> send("dev.off()");
        
        # Make the Location constrained plots
        $R-> send(q`data_for_plot <- cbind(log_otutable_cca_Loc$CCA$wa,log_otutable_cca_Loc$CA$u,mapfile)`);
        #$R-> send("write.table(data_for_plot, file=\"cca_Loc_outdata.txt\", sep=\"\t\", quote=F, col.names=NA)");
        
        $R-> send(q`gp_LocCon_log_Loc <- ggplot() + geom_point(data=data_for_plot, aes(x=CCA1, y=CCA2, color=Tub_Loc, shape=SampTime), size=7)`);
        $R-> send("gp_LocCon_log_Loc_adj <- gp_LocCon_log_Loc + t2 + scale_colour_manual(name=\"Loc\", values=myColors) + labs(title=paste(\"Total constrained: \",CCAperc_log_Loc,\"% Variation\"),x=paste(\"CCA1 (\",CCA1perc_log_Loc,\"% Variation)\"), y =paste(\"CCA2 (\",CCA2perc_log_Loc,\"% Variation)\")) + theme(plot.title = element_text(vjust=3))");
        
        my $outfilename = $resultsdir.$realbase."_LocCon_log_L".$i.".pdf";
        $R-> send("pdf(\"$outfilename\", width=5, height=3.5)");
        $R-> send(q`grid.arrange(gp_LocCon_log_Loc_adj,nrow=1)`);
        $R-> send("dev.off()");
        
        # Make the SampTime constrained plots
        $R-> send(q`data_for_plot <- cbind(log_otutable_cca_SampTime$CCA$wa,log_otutable_cca_SampTime$CA$u,mapfile)`);
        #$R-> send("write.table(data_for_plot, file=\"cca_Loc_outdata.txt\", sep=\"\t\", quote=F, col.names=NA)");
        
        $R-> send(q`gp_SampTimeCon_log_SampTime <- ggplot() + geom_point(data=data_for_plot, aes(x=CCA1, y=CA1, color=Tub_Loc, shape=SampTime), size=7)`);
        $R-> send("gp_SampTimeCon_log_SampTime_adj <- gp_SampTimeCon_log_SampTime + t2 + scale_colour_manual(name=\"Loc\", values=myColors) + labs(title=paste(\"Total constrained: \",CCAperc_log_SampTime,\"% Variation\"),x=paste(\"CCA1 (\",CCA1perc_log_SampTime,\"% Variation)\"), y =paste(\"CA1 (\",CA1perc_log_SampTime,\"% Variation)\")) + theme(plot.title = element_text(vjust=3))");
        
        my $outfilename = $resultsdir.$realbase."_SampTimeCon_log_L".$i.".pdf";
        $R-> send("pdf(\"$outfilename\", width=5, height=3.5)");
        $R-> send(q`grid.arrange(gp_SampTimeCon_log_SampTime_adj,nrow=1)`);
        $R-> send("dev.off()");
        
        $R->stopR();
    }
    
}

close TABLELIST;

exit;

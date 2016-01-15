# Just a simple script.  To use on other directories, change them below.  To use other filters, change them below.
# For values that can be positive and negative (i.e., slopes, filters are for the absolute value)

my $edgedirectory = "./Correlation_Resubmission_Test_NoFilt/All_EdgeTables/";
my $nodedirectory = "./Correlation_Resubmission_Test_NoFilt/All_NodeTables/";
my $tub_supportfilter = "N"; # Comment out line to not filter, otherwise put the value to filter out
my $slopefilter = 0;    # Value of 0 will not filter
my $rsqfilter = 0;      # Value of 0 will not filter
my $rsqadjfilter = 0;   # Value of 0 will not filter
my $pvalfilter = 1;     # Value of 1 will not filter
my $fdrpvalfilter = 1;  # Value of 1 will not filter
my $stmfilter = 0;      # Value of 0 will not filter
my $stmadjfilter = 3;   # Value of 0 will not filter

for($i=4;$i<8;$i++){
    
    my $edgetable = $edgedirectory."Master_EdgeTable".$i.".txt";
    my $nodetable = $nodedirectory."Master_NodeTable".$i.".txt";
    my $edgetable_out = $edgedirectory."Master_EdgeTable".$i."_clean_TubSupp_stmrsqadj3.txt";
    my $nodetable_out = $nodedirectory."Master_NodeTable".$i."_clean_TubSupp_stmrsqadj3.txt";
    
    open(EDGETABLE, $edgetable);
    
    my %edgelinehash;
    while(<EDGETABLE>){
        
        my $currentline = $_;
        chomp $currentline;
        
        my @splitline = split(/\t/, $currentline);
        
        # Skip over correlations based on an organism against itself, these are meaningless - do we need to do this for endo/epi of one with itself?
        if($splitline[0] eq $splitline[1]){next;}                   
        
        
        
        # Sometimes in the tub_support column, a "Y" became a 2 and "N" became a 1 because of R data object issues.  Change this if it is the case
        # In one case it is zero because of some strange entry in the table, I think because of a blank table. Remove this.
        
        if($splitline[2] =~ /[0-9]/){
            if($splitline[2] == 0){next;}
            $splitline[2] =~ s/1/N/g;
            $splitline[2] =~ s/2/Y/g;
            $currentline = join("\t",@splitline);
        }
        
        # Apply any filters that have been set:
        if($splitline[0] =~ /OrgID/){}else{
            if(defined $tub_supportfilter){
                
                if($splitline[2] =~ /$tub_supportfilter/){next;}
            }
            if(abs($splitline[4]) < $slopefilter){next;}
            if(abs($splitline[5]) < $rsqfilter){next;}
            if(abs($splitline[6]) < $rsqadjfilter){next;}
            if($splitline[7] > $pvalfilter){next;}
            if($splitline[8] > $fdrpvalfilter){next;}
            if($splitline[9] < $stmfilter){next;}
            if($splitline[10] < $stmadjfilter){next;}
        }
        
        # Check if the current edge exists in either direction.  If it does, keep the one with the best correlation and ignore the other one.
        
        my $t = 0;
        if(exists $edgelinehash{$splitline[0]}{$splitline[1]}){$t = 1;}
        if(exists $edgelinehash{$splitline[1]}{$splitline[0]}){$t = 2;}
        
        if($t != 0){
            
            if($splitline[0] =~ /OrgID/){                       # If this is just a header line, forget it.
                next;
            }
            
            my $hashkey;
            if($t == 1){
                $hashkey = $edgelinehash{$splitline[0]}{$splitline[1]};
            }
            if($t == 2){
                $hashkey = $edgelinehash{$splitline[1]}{$splitline[0]};
            }
            
            my @splithash = split(/\t/,$hashkey);
            my $stm_rsqadj1 = $splithash[10];
            my $stm_rsqadj2 = $splitline[10];

            if($stm_rsqadj1 > $stm_rsqadj2){
                if($t == 1){
                    $edgelinehash{$splitline[0]}{$splitline[1]} = $hashkey;
                }
                if($t == 2){
                    $edgelinehash{$splitline[1]}{$splitline[0]} = $hashkey;
                }
            }else{
                if($t == 1){
                    $edgelinehash{$splitline[0]}{$splitline[1]} = $currentline;
                }
                if($t == 2){
                    $edgelinehash{$splitline[1]}{$splitline[0]} = $currentline;
                }
            }
        
        }else{
            $edgelinehash{$splitline[0]}{$splitline[1]} = $currentline;
        }
        
    }
    close(EDGETABLE);
    
    open(NODETABLE, $nodetable);
    my %nodelinehash;
    while(<NODETABLE>){
        
        my $currentline = $_;
        chomp $currentline;
        
        my @splitline = split(/\t/, $currentline);
        
        # Now for duplicate nodes, choose BV5 over BV3, Ftrad over FITS2, and Otrad over OITS2 (these choices are not based on anything, except that depth in BV5 is usually more than in BV3)
        if(exists $nodelinehash{$splitline[0]}{$splitline[1]}){
            
            if($splitline[0] =~ /ID/){                       # If this is just a header line, forget it.
                next;
            }
            
            my $hashkey = $nodelinehash{$splitline[0]}{$splitline[1]};
            my @splithash = split(/\t/,$hashkey);
            my $obsin1 = $splithash[3];
            
            my $obsin2 = $splitline[3];
            
            if($obsin1 =~ /^B/){
                
                if($obsin1 =~ /^BV5/){
                    $nodelinehash{$splitline[0]}{$splitline[1]} = $hashkey;
                }else{
                    $nodelinehash{$splitline[0]}{$splitline[1]} = $currentline;
                }
                
            }
            if($obsin1 =~ /^F/){
                
                if($obsin1 =~ /^Ftrad/){
                    $nodelinehash{$splitline[0]}{$splitline[1]} = $hashkey;
                }else{
                    $nodelinehash{$splitline[0]}{$splitline[1]} = $currentline;
                }
                
            }
            if($obsin1 =~ /^O/){
                
                if($obsin1 =~ /^Otrad/){
                    $nodelinehash{$splitline[0]}{$splitline[1]} = $hashkey;
                }else{
                    $nodelinehash{$splitline[0]}{$splitline[1]} = $currentline;
                }
                
            }
        
        
        }else{
            $nodelinehash{$splitline[0]}{$splitline[1]} = $currentline;
        }
        
        
    }
    close(NODETABLE);
    
    # Now we have hashes, let's print them.
    open(EDGETABLEOUT, ">".$edgetable_out);
    open(NODETABLEOUT, ">".$nodetable_out);
    
    my @keys1 = sort { $edgelinehash{$b} cmp $edgelinehash{$a} } keys %edgelinehash;
    foreach (@keys1){
        my $key1 = $_;
        foreach $key2 (keys %{ $edgelinehash{$key1}}) {
            my $printline = "$edgelinehash{$key1}{$key2}\n";
            print EDGETABLEOUT $printline;
        }
    }
    close(EDGETABLEOUT);
    
    my @keys2 = sort { $nodelinehash{$b} cmp $nodelinehash{$a} } keys %nodelinehash;
    foreach (@keys2){
        my $key1 = $_;
        foreach $key2 (keys %{ $nodelinehash{$key1}}) {
            my $printline = "$nodelinehash{$key1}{$key2}\n";
            print NODETABLEOUT $printline;
        }
    }
    close(NODETABLEOUT);
    
}

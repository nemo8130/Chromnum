#!/usr/bin/perl

#==================================
# highest clique finding algorythm
# input : adjasensy matrix representing a directed graph.
#0 1 1 0 0
#1 0 0 1 1
#1 0 0 0 1
#0 1 0 0 1
#0 1 1 1 0
#==================================

#==================================
# read file & Collect data
#==================================


$adjfile = $ARGV[0] || die "Provide the adj file as an argument\n";
chomp $adjfile;

open (INP,"$adjfile");
@adjtest = <INP>;
close INP;

open (OUT,">hclq.report");
open (MAP,">map.chq");

#print scalar(@adjtest),"\n";

@adjm = ();
@nodes = ();

for $a (0..scalar(@adjtest)-1)
{
    @nodes = (@nodes,($a+1));
    @hold = split(/\s+/,$adjtest[$a]);
    $c=0;
    for $b (0..scalar(@hold)-1)
    {
	if ($hold[$b] eq "")
	{
#	    print "I am a null\n";
	}
	elsif ($hold[$b] == 1 || $hold[$b] == 0)
	{
	print $hold[$b]," ";
	$adjm[$a][$c] = $hold[$b];
	$c++;
	}
    }
    print "\n";
}

print "\n";

$tote=0;
$cntS=0;

for $i (0..scalar(@adjm)-1)
{
    for $j (($i+1)..scalar(@adjm)-1)
    {	    
	if ($adjm[$i][$j] == $adjm[$j][$i])
	{
	    $cntS++;
	}
	$tote++;
    }
}

print "$cntS :: $tote\n";

if ($cntS == $tote)
{
    print "The Matrix is symetric :: The graph is directed\n";

    for $i (0..scalar(@adjm)-1)
    {
	$adjm[$i][$i] = 0;              # DIAGONAL ELEMENTS PADDED WITH ZERO
    }

    goto CONTINUE;
}
else
{
    print "The Matrix is not symmetric :: The graph is undirected\n";
    print "Please provide a dirceted graph\n";
    exit;
}

CONTINUE:

@tar = ();
@nei = ();
@cont = ();

print "CONNECTIONS::\n";

for $i (0..scalar(@adjm)-1)
{
    for $j (($i+1)..scalar(@adjm)-1)
    {	    
	print $adjm[$i][$j]," ";
	if ($adjm[$i][$j] == 1)
	{
	    $ti = $i+1;
	    $ni = $j+1;
	    print "$ti -> $ni\n";
	    @tar = (@tar,$ti);
	    @nei = (@nei,$ni);
	    @cont = (@cont,$adjm[$i][$j]);
	}
    }
#   print "\n";
}

#foreach (@tar) {print $_,"\n";}
#print "\n";
#foreach (@nei) {print $_,"\n";}
#print "\n";
#foreach (@nodes) {print $_,"\n";}
#print "\n";

$M = 2;

#==================================
# sort unique nodes
#==================================

%hash1 = ();

foreach $k (@nodes)
{
$hash1{$k}++;
}

@unodes = sort {$a<=>$b} keys %hash1;

foreach (@unodes){print $_,"\n";}

#==================================
# map nodes to index
#==================================

@map = ();
$hash2 = ();

for $i (0..scalar(@unodes)-1)
{
$map[$i] = sprintf("%4d",($i+1)); 
#print $map[$i],"\n";
$hash2{$map[$i]} = $unodes[$i];
}

@mindex = keys %hash2;
@mnodes = values %hash2;

$l = @mindex;
print $l,"\n";

#==================================
# Construct Adjuscency matrix
#==================================

@adj = @adjm;

#=========================================================== check print1
printf  MAP "%3s  %7s     %3s  %7s     %1s\n",'Tar','Tar_map','Nei','Nei_map','adj';

for $i (0..$l-1)
{
	for $j (($i+1)..$l-1)
	{
		if ($adj[$i][$j] == 1)
		{
		printf  MAP "%3d     %4s     %3d     %4s     %1d\n",$map[$i],$unodes[$i],$map[$j],$unodes[$j],$adj[$i][$j];
		}
	}
}

for $i (0..$l-1)
{
#printf "%3d     %7s\n",$mindex[$i],$mnodes[$i];
}
#print "\n\n";
for $i (0..$l-1)
{
#printf "%3d     %7s\n",$map[$i],$unodes[$i];
}

#===============================================

#===============================================
# Map node to degree and neighborhood
#===============================================

@deg = ();
@nh = ();

for $i (0..$l-1)
{
$map[$i];
$deg[$i] = 0;
$c = 0;
	for $j (0..$l-1)
	{
		if ($adj[$i][$j] == 1)
		{
		$deg[$i] = $deg[$i] + $adj[$i][$j];
		$nh[$i][$c] = $map[$j];
		$c++;
		}
	}	
}

#=========================================================== check print2

for $i (0..$l-1)
{
#printf "%3d  %2d\n",$map[$i],$deg[$i];
	for $j (0..$deg[$i]-1)
	{
#	printf "%3d\n",$nh[$i][$j];
	}
#print "\n\n";	
}

#=============================================================
# List triplets
#=============================================================

@clq3 = ();
@clq3map = ();
$c = 0;

for $i (0..$l-1)
{
	for $j (($i+1)..$l-1)
	{
		for $k (($j+1)..$l-1)
		{
			if ($adj[$i][$j] == 1 && $adj[$j][$k] == 1 && $adj[$i][$k] == 1)
			{
			$clq3map[$c][0] = $map[$i];
			$clq3map[$c][1] = $map[$j];
			$clq3map[$c][2] = $map[$k];
			$clq3[$c][0] = $unodes[$i];
			$clq3[$c][1] = $unodes[$j];
			$clq3[$c][2] = $unodes[$k];
			$c++;
			}
		}
	}
}

#=========================================================== check print3

#print $c,"\n";

for $i (0..$c-1)
{
#printf "%7s  %7s  %7s\n",$clq3[$i][0],$clq3[$i][1],$clq3[$i][2];
#printf "%7d  %7d  %7d\n",$clq3map[$i][0],$clq3map[$i][1],$clq3map[$i][2];
}

if ($c >= 1)
{
$M = 3;
}
else
{
print OUT "\n\nNo cliques (of atleast size three) found\n\n";
print  "\n\nNo cliques (of atleast size three) found\n\n";
exit;
}

#============================================================
# The algo
#============================================================

print "=========================================\n";

$N = 3;				# CLIQUE size to start with 

@sortunq = ();

@clqmap = @clq3map;

for $i (0..$c-1)
{
@comp3 = ();
	for $j (0..$N-1)
	{
#	print $clqmap[$i][$j];
	@comp3 = (@comp3,$clqmap[$i][$j]);	
	}
$strcomp = join(' ',@comp3);
print OUT "$strcomp\n";
print "$strcomp\n";
#print "\n";
}

print OUT "\n";
print "\n";

$nnflag = 0;

for $x (3..$l-1)							# highest possible clique size ($l-1) [since it's not a complete graph (tested earlier)
{
$flagN = 0;
@sortunq = ();
$cfreq = 1;
$N = $x;
$nhit = 0;
$nnflag = 0;
	for $i (0..$c-1)
	{
#	print "******************************$cfreq\n";	
	@ano = ();
	$ind = 0;
		for $j (0..$N-1)
		{
		$clqmap[$i][$j];
		@ano = (@ano,$clqmap[$i][$j]);
		}
#	print "clq comp:\n";
#		foreach (@ano){print "\t\t",$_,"\n";}
#	print "\n\n";
	$flag5 = 0;
		for $j (0..scalar(@ano)-1)
		{
		$ind = $ano[$j]-1;
		$deg[$ind];
#		print "clq info: index: $ind    node map: $map[$ind]     deg: $deg[$ind]\n";		
			for $k (0..$deg[$ind]-1)
			{
			$testn = $nh[$ind][$k];			
#			print "test neighbor:            $testn\n";
			$cnt = 0;
#			print "******************",scalar(@ano),"\n";
				for $s (0..scalar(@ano)-1)
				{
				$tit = $ano[$s]-1;
				$tin = $testn-1;
					if ($adj[$tit][$tin] == 1 && $adj[$tin][$tit] == 1)
					{
					$cnt++;
					}
				}
	
				if ($cnt == scalar(@ano))
				{
				$flag5 = 1;
#				print "nclq found for: \n";
					foreach (@ano)
					{
#					print $_,"\n";
					}
#				print "\nnew mem: $testn\n";
				$nhit++;							
#				print "new feed:                                                        $testn\n";
				$clqmap[$i][scalar(@ano)] = $testn;			
#				print "New clique found for: old clique $cfreq\n";
#				print "\n\n";
				$flagN = 1;
				$j = scalar(@ano)-1;
				$M = $N+1;
					for $a (0..$nhit-1)
					{
						for $b (0..$M-1)
						{
						@store = ();
							for $j (0..scalar(@ano))
							{
#								if ($clqmap[$i][$j] != 0)
#								{
								$hgclq[$a][$b] = $clqmap[$i][$j];
#								print "****************CHQ******************* $hgclq[$a][$b]\n";
								@store = (@store,$hgclq[$a][$b]);
								$nnflag = 1;
#								}
							}
#							if ($nnflag == 1)
#							{
							@astore = sort {$a<=>$b} @store;
							$strn = join('',@astore);
#							}							
						}
					}
#					if ($nnflag == 1)
#					{
					@sortunq = (@sortunq,$strn);
#					}
				}
			}
		}
	$cfreq++;
	}

#	print "===================================~~~~~~~~~~~~~~~~\n";
#	foreach (@sortunq){print $_,"\n";}
#	print "===================================~~~~~~~~~~~~~~~~\n";

#print "===============================flag = ",$flagN,"\n";
	if ($flagN == 0)
	{
	print OUT "Nodes involved in the largest clique :\n";
	print "Nodes involved in the largest clique :\n";
		for $u (0..$c-1)
		{
			for $v (0..$x-1)
			{
#			print $clqmap[$u][$v],"  ";
			print OUT $clqmap[$u][$v],"  ";
			print  $clqmap[$u][$v],"  ";
			}
		print "\n";
		print OUT "\n";
		}
	last;
	}
	elsif ($flagN == 1)
	{
	%shash = ();

		foreach $z (@sortunq)
		{
		$shash{$z}++;
		}
		@ncomp = ();
		@ncomp = keys %shash;
	
		$lnn = @ncomp;
#		print "                         $lnn\n";
	
		foreach (@ncomp)
		{
#		print "here it is : $_\n";
		print OUT "$_\n";
		print  "$_\n";
		}
	print OUT "\n\n";
	print  "\n\n";
	}
	
	
	
@nnc = ();

$col = 0;
$row = 0;	
	
	for $u (0..$lnn-1)
	{
	$col = 0;
		for ($v=0;$v<length($ncomp[$u]);$v+=4)
		{
		$nnc[$row][$col] = sprintf("%4d",(int(substr($ncomp[$u],$v,4))));
#		print "######################",$nnc[$row][$col],"\n";
		$col++;
		}
	$row++;
	}
	
#print "$x   $M    $lnn    $c","   ",length($ncomp[$u]),"\n";
	
@clqmap = ();
	
	for $u (0..$lnn-1)
	{
		for $v (0..$M-1)
		{	
#		print "-------------------$nnc[$u][$v]\n";
		$clqmap[$u][$v] = $nnc[$u][$v];
#		print "=============================\n";
#		print "$nnc[$u][$v]       <-               $clqmap[$u][$v]\n";
#		print "=============================\n";
		}
	}	
$c = $lnn;
$N = $M;	
}	



if ($flagN == 1)
{
print OUT "Nodes involved in the largest clique :\n";
	for $u (0..$c-1)
	{
		for $v (0..$M-1)
		{
#		print $clqmap[$u][$v],"  ";
		print OUT $clqmap[$u][$v],"  ";
		print  $clqmap[$u][$v],"  ";
		}
	print "\n";
	print OUT "\n";
	}
}
print OUT "\n\n";
print OUT "Size of the largest clique : $M\n";
print OUT "Number of cliques of the highest order: $c\n";

print "\n\n";
print "Size of the largest clique : $M\n";
print "Number of cliques of the highest order: $c\n";
			

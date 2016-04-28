#!/tools/bin/perl
use strict;
use lib '/home/people/rachita/perl5/lib/perl5/ia64-linux-thread-multi-ld/auto';
use lib '/home/people/rachita/perl5/lib';
use DBI;
use List::MoreUtils qw/ uniq /;
#use Tree ();
use Data::Dumper;

my $results = $ARGV[0];	## for GO lists with coln in order ID gene list name number p-val tab seperated
#my $number = $ARGV[1];	## top hits you want
#my $output = $ARGV[2]; ##OUTPUT FILE

 my $dbh = DBI->connect("DBI:mysql:database=rachita_private;host=mysql.cbs.dtu.dk","rachita","PfbMOOay", {'RaiseError' => 1}) or die "Connection Error: $DBI::errstr\n";
 
 open(INGO,'<',$results) or die;
 my %GO_results=();
 
 while(defined(my $line =<INGO>))
 {
 	chomp($line);
	my @row=split("\t",$line);
 	$GO_results{$row[0]}=[@row];
 }
 
my @class_set=keys %GO_results;
my %child_parent=();
my %parents=();

 foreach my $key (keys %GO_results)
 {
 	print $key,"\n";
	
	my $sql = "SELECT DISTINCT
        term.acc,term.name,
        ancestor.acc, ancestor.name,
        graph_path.distance
	FROM 
  	go_termdb.term
  	INNER JOIN go_termdb.graph_path ON (term.id=graph_path.term2_id)
 	INNER JOIN go_termdb.term AS ancestor ON (ancestor.id=graph_path.term1_id)
 	WHERE term.acc = ?
	and graph_path.distance <> 0
	and ancestor.acc in (" . join(",", map { $dbh->quote($_) } @class_set). ")
	order by graph_path.distance;";
	#print "$sql\n";
 	my $sth = $dbh->prepare($sql);
	$sth->bind_param( 1, $key );
	#print "$sth\n";
 	$sth->execute or die "SQL Error: $DBI::errstr\n";
	#print "$sth\n",$sth->rows,"\n";
	#print $sth->fetchrow_array;
	
 	while (my @row = $sth->fetchrow_array) {
		print @row,"\n";
 		my $child=$row[0];
		my $parent=$row[2];
		my $distance=$row[4];
		my $par_name=$row[3];
		
		my $sql1 = "SELECT DISTINCT
   		term.*,p.distance
 		FROM
   		go_termdb.term
   		INNER JOIN go_termdb.graph_path AS p ON (p.term2_id=term.id)
   		INNER JOIN go_termdb.term AS root ON (p.term1_id=root.id)
 		WHERE
   		root.is_root=1 and term.acc=?";
	print "$sql\n";
 		my $sth1 = $dbh->prepare($sql1);
		$sth1->bind_param( 1, $parent );
	#print "$sth\n";
 		$sth1->execute or die "SQL Error: $DBI::errstr\n";
		while (my @row1 = $sth1->fetchrow_array) {
			print "Parent:", $row1[1],"\t",$row1[4],"\n";
		}
		if (exists($parents{$parent})){
			$parents{$parent}=$parents{$parent}+1;
		}
		else
		{
			$parents{$parent}=1;
		}
		my @children=@{$child_parent{$parent}{'child'}} if exists($child_parent{$parent}{'child'});
		push(@children,$child);
		$child_parent{$parent}{'child'}=[@children];
		$child_parent{$parent}{'distance'}=$distance;
 	} 

 
 }
close(INGO);

#my $max_class=(sort {$parents{$b} <=> $parents{$a}} keys %parents)[1];
#print $max_class,"\n";

open(OUT,'>',"CYTOSCAPE_file.sif") or die;

for my $key(sort {$parents{$b} <=> $parents{$a}} keys %parents)
{
	print $key,":",$parents{$key},"-", uniq(@{$child_parent{$key}{'child'}}),"\n";
	#print $key," = ";
	foreach my $ch (uniq(@{$child_parent{$key}{'child'}}))
	{
		print OUT $key," = ",$ch,"\n";
	}
}

close(OUT);

#open(CP,'<',"CYTOSCAPE_file") or die;
#open(OX,'>',"par_child") or die;
#
#my %nodes;
#my @children; 
#while (<CP>) {
#   my ($c, $p) = split /=/;
#  push @children,$c unless $nodes{$c};
#    $nodes{$_}{name}||=$_
#        for $c,$p;
#    $nodes{$p}{kids}{$c}=$nodes{$c};
#}
#delete $nodes{$_}
#    for @children;
#my @roots=keys %nodes;
#
#print OX "\nHERE $_",print_inorder($nodes{$_}) for @roots;
#
#sub print_inorder {
#   my $node= shift;
#   print OX "THESE", $node->{name}, " ";
#   print_inorder($_)
#      for sort { $a->{name} cmp $b->{name} } 
#          values %{ $node->{kids} || {} };
#}
#
##and ancestor.acc in ('')
##and ancestor.acc in (?)
##$sth->bind_param( 2, $all_classes);

# script accompanying the paper 
# polydb: A Database for Polytopes and related Objects
# 
# functions:
# identify_smooth_fano_in_polydb
# all_free_sums_in_dim
# skew_sum_with_simplex
# skew_simplex_k_sums_in_dim
# all_skew_simplex_sums_in_dim
# all_skew_bipyramids_in_dim


use application "polytope";

# checks for the name/id of a smooth Fano polytope
# in the database polyDB
# input: $p Polytope
# return: id in polydb
sub identify_smooth_fano_in_polydb {
	my $p = shift;	
	my $d = $p->DIM;
	my $nlp = new Int($p->N_LATTICE_POINTS);
	my $parray=db_query({"DIM"=>$d, "N_VERTICES"=>$p->N_VERTICES, "N_FACETS"=>$p->N_FACETS, 
		"N_LATTICE_POINTS"=>$nlp, }, db=>"LatticePolytopes", collection=>"SmoothReflexive");
	foreach my $c ( @$parray ) {
		if ( lattice_isomorphic_smooth_polytopes($c,$p) ) { return $c->name; }
	}
	die "polytope not found\n";
}

# determines ids of all free sums among the smooth Fano polytopes of a given dimension
# input $d: dimension
#       %options: verbose: print progress info, splitinfo: additionally return all possible decompositions for an id
# return list of ids or map assinging each id its possible decompositions
sub all_free_sums_in_dim { 
	if ( @_ == 0 ) {
		print "usage: all_free_sums_in_dim(dimension, options) where allowed options are verbose, splitinfo, start, end, skip, amount given as hash, e.g. verbose=>1, skip and amount are considered wrt to the second factor, i.e. for dimension 6, start=end=2, amount=10, skip=50 we get all products of a 2D with a 4D polytope where the index of the 4D polytope is between 50 and 59";
		return;
	}

	my ($d,%options) = @_;
	my $list;
	if ( $options{"splitinfo"} ) {
		$list = new Map<String,Set<Pair<String,String> > >;
	} else {
		$list = new Set<String>;
	}
	
	my $start = 1;
	my $end = $d/2;
	if ( $options{"start"} ) {
		$start = $options{"start"};
	}
	if ( $options{"end"} ) {
		$end = $options{"end"};
	}
	
	my $cur_options1 = { db=>"LatticePolytopes", collection=>"SmoothReflexive" };
	if ( $options{"skip"} ) {
		$$cur_options1{"skip"} = $options{'skip'};
	}
	if ( $options{"amount"} ) {
		$$cur_options1{"limit"} = $options{'amount'};
	}
	my $cur_options2 = { db=>"LatticePolytopes", collection=>"SmoothReflexive" };
	
	if ( $options{"verbose"} ) {
		print "checking for splits of a $d-dimensional polytope\n";
	}
	foreach my $n ($start..$end) {
		if ( $options{"verbose"} ) {
			print "checking for splits into a $n-dimensional and a ", $d-$n, "-dimensional polytope\n";
		}
		my $cur1=db_cursor({"DIM"=>$d-$n}, $cur_options1 );
		while ( !$cur1->at_end() ) {
			my $c1 = $cur1->next();
			my $cur2=db_cursor({"DIM"=>$n}, $cur_options2 );
			while ( !$cur2->at_end() ) {
				my $c2 = $cur2->next();
				next if ( $n == $d-$n && $c2->name lt $c1->name );
				my $name = identify_smooth_fano_in_polydb(product($c1,$c2));
				if ( $options{"verbose"} ) {
					print "split found: $name splits into ", $c1->name, " and ", $c2->name, "\n";
				}
				if ( $options{"splitinfo"} ) {
					my $split = new Pair<String,String>($c1->name,$c2->name);
					$list->{$name} += $split;
				} else {
					$list += $name;
				}
			}
		}
	}
	return $list;
}


# returns the dual of the skew sum of the polytope with a simplex, where one vertex is shifted
# input: $q: smooth Fano polytope
#        $k: dimension of the simplex
#        $n: the shifted vertex
sub skew_sum_with_simplex {
	my ($q,$k,$n) = @_;
	my $m = new Matrix(primitive($q->FACETS));
	$m |= (zero_matrix($m->rows,$k));
	foreach (1..$k) {
		$m /= -(unit_vector($m->cols,$m->cols-$_));
		$m->elem($m->rows-1,0) = 1;
	}
	$m /= $n;
	return new Polytope(INEQUALITIES=>$m);
}

# determines ids of all generalized smooth simplex sums in dimension d with a simplex of dimension k
# input $d: dimenison 
#       $k dimension of simplex
#       %options: 
#          verbose: print progress
#          splitinfo: additionally return all decompositions
#          skip: omit the first skip smooth Fanos of dimension d-k from the db (to recover an interupted computation)
#          amount: only process amount smooth Fanos in dimension d-k from the db (to recover an interupted computation)
sub skew_simplex_k_sums_in_dim {

	if ( @_ == 0 ) {
		print "usage: skew_simplex_k_sums_in_dim(dimension, simplex_dimension, options) where allowed options are skip, amount, verbose, splitinfo, given as hash, e.g. skip=>1000";
		return
	}
	
	my ($d,$k,%options) = @_;
	die "simplex dimension too big" if  $d < $k;
	die "dimension too big" if $d > 9;

	my $sp;
	if ( $options{"splitinfo"} ) {
		$sp = new Map<String,Set<Pair<String,Vector<Integer>>>>;
	} else {
		$sp = new Set<String>;
	}


	my $simplex = new Matrix<Integer>($k+1,$d+1);
	$simplex->col(0) = ones_vector<Integer>($k+1);
	foreach (1..$k) {
		$simplex->elem($_,$d-$k+$_) = -1;
		$simplex->elem(0,$d-$k+$_) = 1;
	}
	
	if ( $options{"verbose"} ) {
		print "checking simplex sums in dimension $d with simplex dimension $k ";
		if ( $options{"skip"} ) {
			print "skipping first $options{'skip'} ";
		}
		if ( $options{"amount"} ) {
			print "limited to $options{'amount'} polytopes";
		}
		print "\n";
	}
	
	if ( $d == $k ) {
		my $r = new Polytope(FACETS=>$simplex);
		my $rname = identify_smooth_fano_in_polydb($r);
		my $n = $simplex->row(0);
		if ( $options{"verbose"} ) {
			print "found a new simplex sum: ", $rname, " is the sum of the origin with a $k-dimensional simplex, one vertex shifted into ", $n, "\n";
		}
		if ( $options{"splitinfo"} ) {
			$sp->{$rname} += new Pair<String,Vector<Integer>>($rname,$n);
		} else {
			$sp += $rname;
		}
		return $sp;
	}
		
	my $eq = new Matrix<Integer>($k,$d+1);
	foreach (0..$k-1) {
		$eq->elem($_,0) = 1;
		$eq->elem($_,$d-$k+$_+1) = -1;
	}
	my $cur_options = { db=>"LatticePolytopes", collection=>"SmoothReflexive" };
	if ( $options{"skip"} ) {
		$$cur_options{"skip"} = $options{'skip'};
	}
	if ( $options{"amount"} ) {
		$$cur_options{"limit"} = $options{'amount'};
	}

	my $cur=db_cursor({"DIM"=>$d-$k}, $cur_options);
	while ( !$cur->at_end() ) {
		my $p = $cur->next;
		my $m = new Matrix<Integer>(primitive($p->FACETS));
		$m |= zero_matrix<Integer>($m->rows,$k);
		my $s = new Polytope(POINTS=>$simplex/$m);
		my $mm = $s->FACETS->minor(~$s->FACETS_THRU_VERTICES->[0],All);
		my $q=new Polytope(INEQUALITIES=>$mm, EQUATIONS=>$eq);
		my $lp = $q->INTERIOR_LATTICE_POINTS;
		
		if ( $options{"verbose"} ) {
			print $p->name, "---------------- \n";
			print "checking the points:\n", $lp, "\n";
		}
		foreach my $n (@$lp) {
			my $r = skew_sum_with_simplex($p,$k,$n);
			my $rname = identify_smooth_fano_in_polydb($r);
			if ( $options{"verbose"} ) {
				print "found a new simplex sum: ", $rname, " is the sum of ", $p->name, " with a $k-dimensional simplex, one vertex shifted into ", $n, "\n";
			}
			if ( $options{"splitinfo"} ) {
				$sp->{$rname} += new Pair<String,Vector<Integer>>($p->name,$n);
			} else {
				$sp += $rname;
			}
		}
	}
	return $sp;
}


# returns ids of all generalized simplex sums in a given dimension
# input: $d: dimension
#        %options: verbose: print progress information, splitinfo: also return all possible decompositions
# returns list of ids, or Map that assigns an id its possible decompositions, if splitinfo is set
sub all_skew_simplex_sums_in_dim {
	if ( @_ == 0 ) {
		print "usage: all_skew_simplex_sums_in_dim(dimension, options) where allowed options are verbose, splitinfo, given as hash, e.g. skip=>1000";
		return
	}
	my ($d,%options) = @_;
	die "dimension too big" if $d > 9;

	my $sp;
	if ( $options{"splitinfo"} ) {
		$sp = new Map<String,Set<Pair<String,Vector<Integer>>>>;
	} else {
		$sp = new Set<String>;
	}
	
	foreach ( 1..$d ) {
		my $result = skew_simplex_k_sums_in_dim($d,$_,%options);
		if ( $options{"splitinfo"} ) {
			foreach my $key (keys %$result) {
				$sp->{$key} += $result->{$key};
			}
		} else {
			$sp += $result;
		}
	}
	
	return $sp;
}


# determines all skew bipyramids in a given dimension 
# according to the definition in Assarf et al., Smooth Fano polytopes with many vertices
# input: $d dimension
# options:
#	verbose: default 0, if set to 1 then print some information on the progress 
#	splitinfo: default 0, in this case returns a Set with names/ids of polytopes that are a skew bipyramid
# 		if set to 1, then the function returns a Map assigning each id of a polytope 
# 		that is a skew bipyramid a list of all potential based of the bipyramid and the vertices above which one apex is chosen
#	skip n: skip first n polytopes in idm d-1
#	amount n: consider only n polytopes in dimension d-1
sub all_skew_bipyramids_in_dim { 
	
	if ( @_ == 0 ) {
		print "usage: all_skew_bipyramids_in_dim(dimension, options) where allowed options are verbose, splitinfo, skip, amount, given as hash, e.g. verbose=>1";
		return
	}
		
	my ($d,%options) = @_;
	my $list;
	if ( $options{"splitinfo"} ) {
		$list = new Map<String,Set<Pair<String,Int> > >;
	} else {
		$list = new Set<String>;
	}

	if ( $options{"verbose"} ) {
		print "checking for skew bipyramids among the $d-dimensional smooth Fano polytopes ";
		if ( $options{"skip"} ) {
			print "skipping first $options{'skip'} ";
		}
		if ( $options{"amount"} ) {
			print "limited to $options{'amount'} polytopes";
		}
		print "\n";
	}
	
	my $cur_options = { db=>"LatticePolytopes", collection=>"SmoothReflexive" };
	if ( $options{"skip"} ) {
		$$cur_options{"skip"} = $options{'skip'};
	}
	if ( $options{"amount"} ) {
		$$cur_options{"limit"} = $options{'amount'};
	}
	
	my $cur=db_cursor({"DIM"=>$d-1}, $cur_options);
	my $bottom = dense(unit_vector<Integer>($d+1,0));
	$bottom->[$d] = -1;
	while ( !$cur->at_end() ) {
		my $c = $cur->next();
		my $fac = primitive($c->FACETS);
		for my $i (0..$c->N_FACETS-1) {
			my $top = new Vector<Integer>($fac->row($i)|1);
			my $p = new Polytope(FACETS=>($fac|zero_vector<Integer>($c->N_FACETS))/$top/$bottom);
			my $name = identify_smooth_fano_in_polydb($p);
			if ( $options{"verbose"} ) {
				print "skew bipyramid found: $name is a skew bipyramid over ", $c->name, " and vertex ", $i, "\n";
			}
			if ( $options{"splitinfo"} ) {
				$list->{$name} += new Pair<String,Int>($c->name,$i);
			} else {
				$list += $name;
			}
		}
	}
	return $list;
}


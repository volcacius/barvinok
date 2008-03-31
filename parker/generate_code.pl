#!/usr/local/bin/perl

# Erin Parker (parker@cs.unc.edu), March 2004

# generate_code.pl  -- Takes as input a Presburger formula and gives as output source code 
#                      for building a DFA representation of the formula and subsequently 
#                      counting the number of accepting paths in the DFA to reveal the 
#                      number of solutions to the original formula.  (For details on 
#                      this method see Parker & Chatterjee, Compiler Construction 2004.)

# The input formula is assumed to be written in the style of the Omega Library.  The 
# output source code makes calls to Constantinos Bartzis' DFA-construction code and the
# MONA tool's DFA Library.  This code is certainly not guaranteed to be bug-free.  It really 
# is designed to process the sort of representations I typically see when using the Omega 
# Library for the expression and manipulation of Presburger formulas of interest to me. 
#                                                                                   ---Erin

# FORMULA:     "{[TUPLE] : CLAUSE union 
#               CLAUSE union 
#               ... union 
#               CLAUSE}"
# CLAUSE:      "CONSTRAINT && CONSTRAINT && ... && CONSTRAINT"
# CONSTRAINT:  "coeff*var +/- ... (+/- const) eq/ineq coeff*var +/- ... (+/- constant)"   

# NumFreeVars -- number of free variables, given on command line
# FreeVar -- array of free variable names, "In_1 In_2 . . . In_{NumFreeVars}"
# ClauseMinIndex, ClauseMaxIndex -- CLAUSE i owns CONSTRAINTs 
#   ClauseMinIndex to ClauseMaxIndex
# NumStrideVars -- for CLAUSE i, NumStrideVars gives number of stride variables
# ConstraintMinIndex, ConstraintMaxIndex -- CONSTRAINT i owns variables 
#   (from array Var) ConstraintMinIndex to ConstrainMaxIndex with 
#    coefficients (from array Coeff)
# Constant -- CONSTRAINT i owns constant i
# Kind -- CONSTRAINT i is of type $Eq or $Ineq   
# MaxConst -- MaxConst i is the largest coeff or constant in CONSTRAINT i 


$False = 0;
$True = 1;

$Eq = 0;
$Ineq = $LTE = 1;
$Less = 2; 


# SUBROUTINE: StoreConstraint($Kind, $LHS, $RHS, $NumVars)  
sub StoreConstraint {

  my($Kind);
  $Kind = $_[0];
  my ($LHS);
  $LHS = $_[1];
  my ($RHS);
  $RHS = $_[2];
  my ($NumVars);
  $NumVars = $_[3];

  $ConstraintMinIndex[$ConstraintNum] = $VarNum;
  $MaxConst[$ConstraintNum] = 0;

  $LHS =~ s/\+/\$\+/g;
  $LHS =~ s/\-/\$\-/g;
  $RHS =~ s/\+/\$\+/g;
  $RHS =~ s/\-/\$\-/g;

  if($Kind == $LTE) {
    $Constant[$ConstraintNum] = -1;
    $Kind[$ConstraintNum] = $Ineq;
  }
  else {
    $Constant[$ConstraintNum] = 0;
    if($Kind == $Eq) {
      $Kind[$ConstraintNum] = $Eq;
    }
    else {
      $Kind[$ConstraintNum] = $Ineq;
    }
  }  

  @OperandList = split(/\$/, $LHS);  
  foreach $Operand (@OperandList) { 
    if($Operand ne "") {  
      $Operand =~ s/\+//; 
      $Found = $False;
      for($x=0; $x<$NumVars; $x++) { 
        if($Operand =~ /$FreeVar[$x]/ && $Operand !~ /$FreeVar[$x]\S/ &&
	   ($Operand !~ /\D$FreeVar[$x]/ || $Operand =~ /-$FreeVar[$x]/)) {
	  $Operand =~ s/$FreeVar[$x]//;
          $Found = $True;
	  $Index = $x;
	  $x = $NumVars;
        }
      } 
      if($Operand > $MaxConst[$ConstraintNum]) {
	$MaxConst[$ConstraintNum] = $Operand;
      }
      if(!$Found) {
        $Constant[$ConstraintNum] += $Operand;
      }
      else {
        if($Operand eq "") {
	   $Operand = 1;
        }
	elsif($Operand eq "-") {
	    $Operand = -1;
	}
        $Var[$VarNum] = $Index;
        $Coeff[$VarNum++] = $Operand;
      }
    }
  }

  @OperandList = split(/\$/, $RHS);
  foreach $Operand (@OperandList) {
    if($Operand ne "") {   
      $Operand =~ s/\+//;
      $Found = $False;
      for($x=0; $x<$NumVars; $x++) {
	if($Operand =~ $FreeVar[$x] && $Operand !~ /$FreeVar[$x]\S/ &&
	   ($Operand !~ /\D$FreeVar[$x]/ || $Operand =~ /-$FreeVar[$x]/)) {
	  $Operand =~ s/$FreeVar[$x]//;
	  $Found = $True;
	  $Index = $x;
	  $x = $NumVars;
	}
      }
      if($Operand > $MaxConst[$ConstraintNum]) {
	$MaxConst[$ConstraintNum] = $Operand;
      }
      if(!$Found) {
	if($Operand != 0) {
	  $Operand *= -1;
        }
	$Constant[$ConstraintNum] += $Operand;
      }
      else {
	if($Operand eq "") {
	    $Operand = 1;
	}
	elsif($Operand eq "-") {
	    $Operand = -1;
	}
        $Var[$VarNum] = $Index;
        $Coeff[$VarNum++] = -1*$Operand;
      }
    }
  }
  $ConstraintMaxIndex[$ConstraintNum++] = $VarNum-1;
}


# SUBROUTINE: ParseAndStore()
sub ParseAndStore {

  # check for correct command line usage
  if(@ARGV < 1) {
    print "USAGE: num_free_vars \n\n";
    exit;
  }

  # create free variable names
  $NumFreeVars = $ARGV[0];  
  for($z=0; $z<$NumFreeVars; $z++) {
    $FreeVar[$z] = "In_".($z+1);
  }

  print "/* This code is automatically generated to build a DFA\n"; 
  print "   representing integer solutions to the following formula.\n\n";

  $ClauseNum = 0;
  $ConstraintNum = 0;
  $MaxConstraintNum = 0;
  $VarNum = 0;

  $NumLines = 0;

  while($Line = <STDIN>) {

    if($Line !~ /\#/ && $Line ne "") {
      print "   ".$Line;
    }      

    # ignore whitespace, braces and "union" 
    $Line =~ s/[\s|\}|\{]//g;
    $Line =~ s/union//g;

    # avoid Omega comments and blank lines  
    if($Line !~ /\#/ && $Line ne "") {

      $Skip = 0;  
       for($q=0; $Skip==0 && $q<$NumLines; $q++) {
        if($Line eq $LineList[$q]) {
	  $Skip = 1;
        }
      }
      if($Skip==0) {
       $LineList[$NumLines++] = $Line;


      $ClauseMinIndex[$ClauseNum] = $ConstraintNum;


      # split TUPLE from the rest of the line
      ($Tuple, $Rest) = split("]:", $Line);
      $Tuple =~ s/\[//g;  
      @TupleVar = split(",", $Tuple);

      # create any necessary constraints from list of tuple variables
      for($i=0; $i<=$#TupleVar; $i++) {
    
        if($TupleVar[$i] ne $FreeVar[$i]) {
          StoreConstraint($eq, $FreeVar[$i], $TupleVar[$i], $NumFreeVars);
        }
      }

      # split the rest of the line into CONSTRAINTs
      @ConstraintList = split("&&", $Rest); 
  
      # handle "Exists(..." in input formula
      $NumStrideVars[$ClauseNum] = 0;
      if($ConstraintList[0] =~ /Exists/) {
	$ConstraintList[0] =~ s/Exists\(//g;
	($ExistsPart, $ConstraintList[0]) = split(":", $ConstraintList[0]);
	@StrideVar = split(",", $ExistsPart);

        for($e=0; $e<=$#StrideVar; $e++, $NumStrideVars[$ClauseNum]++) {
          $FreeVar[$NumFreeVars+$e] = $StrideVar[$e];
        }    

	$ConstraintList[-1] =~ s/\)//g;    
      }

      $TotalVars = $NumFreeVars + $NumStrideVars[$ClauseNum];

      # process each CONSTRAINT
      foreach $Constraint (@ConstraintList) { 

        # handle '<=' in CONSTRAINT
        if($Constraint =~ /<=/) {
	  $Constraint =~ s/<=/\$/g;
	  # Note: This does not handle every possible combination of 
	  #       '<=' and '<' in constraint.
	  $Lspecial = $False;
	  $Rspecial = $False;
	  ($Lside, $Rside) = split(/\$/, $Constraint);
	  if ($Lside =~ /</) {
	    $Lspecial = $True;
	    $Constraint =~ s/</\$/g;
	  }
	  elsif ($Rside =~ /</) {
	    $Rspecial = $True;
	    $Constraint =~ s/</\$/g;
	  }

          @SubConstraint = split(/\$/, $Constraint);
          if($#SubConstraint > 3) {
	    print "Too many \"<=\" and \"<\" in a single constraint\n\n";
            exit;
          }
          @Part0 = split(/,/, $SubConstraint[0]);
          @Part1 = split(/,/, $SubConstraint[1]);
          @Part2 = split(/,/, $SubConstraint[2]);
          @Part3 = split(/,/, $SubConstraint[3]);
      
          foreach $LTE_Lexpr (@Part0) {
            foreach $LTE_Rexpr (@Part1) {
	      if ($Lspecial) {
                StoreConstraint($Less, $LTE_Lexpr, $LTE_Rexpr, $TotalVars);
	      }
	      else {
		StoreConstraint($LTE, $LTE_Lexpr, $LTE_Rexpr, $TotalVars);
	      }
	    }
          }   

          foreach $LTE_Lexpr (@Part1) {
            foreach $LTE_Rexpr (@Part2) {
	      if ($Rspecial) {
                StoreConstraint($Less, $LTE_Lexpr, $LTE_Rexpr, $TotalVars);
	      }
	      else {
                StoreConstraint($LTE, $LTE_Lexpr, $LTE_Rexpr, $TotalVars);
	      }
            }
          }     

          foreach $LTE_Lexpr (@Part2) {
            foreach $LTE_Rexpr (@Part3) {
              StoreConstraint($LTE, $LTE_Lexpr, $LTE_Rexpr, $TotalVars);
            }
          } 
        }  

        # handle '=' in CONSTRAINT
        elsif($Constraint =~ /=/) {
          @SubConstraint = split(/=/, $Constraint);
        
	  @Part0 = split(/,/, $SubConstraint[0]);
	  @Part1 = split(/,/, $SubConstraint[1]);
      
	  foreach $EQ_Lexpr (@Part0) {
	    foreach $EQ_Rexpr (@Part1) {
	      StoreConstraint($Eq, $EQ_Lexpr, $EQ_Rexpr, $TotalVars);
            }
          }   
        }

        # handle '<' in Constraint
        elsif($Constraint =~ /</) {
          @SubConstraint = split(/</, $Constraint);
        
	  @Part0 = split(/,/, $SubConstraint[0]);
	  @Part1 = split(/,/, $SubConstraint[1]);
	  @Part2 = split(/,/, $SubConstraint[2]);
 
          if($#SubConstraint > 2) {
	    print "Too many \"<\" in a single constraint\n\n";
            exit;
          }
     
          foreach $LESS_Lexpr (@Part0) {
            foreach $LESS_Rexpr (@Part1) {
              StoreConstraint($Less, $LESS_Lexpr, $LESS_Rexpr, $TotalVars);
            }
          }   

          foreach $LESS_Lexpr (@subLESSs_1) {
            foreach $LESS_Rexpr (@subLESSs_2) {
              StoreConstraint($Less, $LESS_Lexpr, $LESS_Rexpr, $TotalVars);
            }
          }   
        }
      }
      if(2*($ConstraintNum - $ClauseMinIndex[$ClauseNum] - 1) - 1 > $MaxConstraintNum) {
	$MaxConstraintNum = 2*($ConstraintNum - $ClauseMinIndex[$ClauseNum] - 1) - 1;
      }
      $ClauseMaxIndex[$ClauseNum++] = $ConstraintNum-1; 
    }
   }
  }
  print "*/\n\n";
}


# SUBROUTINE GenerateCode()
sub GenerateCode {

  print "#include \"count_solutions.h\"\n\n";

  print "int main(void) {\n\n";
  print "  int coeffs[".($NumFreeVars+3)."];\n";       
  print "  int indices[".($NumFreeVars+3)."];\n";
  print "  DFA* dfa1[".($MaxConstraintNum+4)."];\n";
  print "  DFA* dfa2[".(2*$ClauseNum-1)."];\n\n";

  for($i=0, $n=0; $i<$ClauseNum; $i++, $n++) {
    for($j=$ClauseMinIndex[$i], $m=0; $j<=$ClauseMaxIndex[$i]; $j++, $m++) {
      for($k=$ConstraintMinIndex[$j], $l=0; $k<=$ConstraintMaxIndex[$j]; $k++, $l++) {
        print "  coeffs[".$l."] = ".$Coeff[$k].";\n";
	print "  indices[".$l."] = ".$Var[$k].";\n";
      }
      print "  dfa1[".$m."] = build_DFA_";
      if($Kind[$j]) {
	print "in";
      }
      print "eq(".$l.", coeffs, ".$Constant[$j].", indices);\n";
      if($m>0) {
        print "  dfa1[".(++$m)."] = dfaMinimize(dfaProduct(dfa1[".($m-2)."], dfa1[".($m-1)."], dfaAND));\n";
	print "  dfaFree(dfa1[".($m-2)."]);\n  dfaFree(dfa1[".($m-1)."]);\n";
      }
      print "\n";
    }

    for ($f=0; $f<$NumStrideVars[$i]; $f++) {
      if ($f == $NumStrideVars[$i]-1) {
	  print "  dfa2[".($n)."] = dfaMinimize(dfaProject(dfa1[".($m-1)."], ".($NumFreeVars+$f)."));\n";
	  print "  dfaFree(dfa1[".($m-1)."]);\n";
      }
      else {
	print "  dfa1[".($m++)."] = dfaMinimize(dfaProject(dfa1[".($m-2)."], ".($NumFreeVars+$f)."));\n";
	print "  dfaFree(dfa1[".($m-2)."]);\n";
      }
    }

    if ($NumStrideVars[$i] == 0) {
      print "  dfa2[".($n)."] = dfaCopy(dfa1[".($m-1)."]);\n";
      print "  dfaFree(dfa1[".($m-1)."]);\n";
    }

    if($n>0) {
      print "  dfa2[".(++$n)."] = dfaMinimize(dfaProduct(dfa2[".($n-2)."], dfa2[".($n-1)."], dfaOR));\n";
      print "  dfaFree(dfa2[".($n-2)."]);\n  dfaFree(dfa2[".($n-1)."]);\n";
    }
    print "\n\n";
  }

  print "  count_accepting_paths(dfa2[".(2*$ClauseNum-2)."], dfa2[".(2*$ClauseNum-2)."]->ns, ".$NumFreeVars.");\n";
  print "  dfaFree(dfa2[".(2*$ClauseNum-2)."]);\n";
  print "}\n\n";
}


# MAIN ROUNTINE

ParseAndStore();
GenerateCode();









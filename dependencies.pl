#! /usr/bin/perl
#

my $computer = "office";
my %params;

my $HOME = `cd; pwd`; chomp $HOME;

  @{$params{"defs"}} = ( "_G_HAVE_BOOL", "ANSI_HEADERS", "HAVE_SSTREAM" ); 

@{$params{"Includes"}} = ( "." );

if ( $computer =~ /office/ )
{
  @{$params{"make include"}} = ( "/usr/local/include", 
                                 "$HOME/usr/include",
                                 "/usr/local/include/newmat",
                                 "/usr/local/include/opt++",
                                 "/usr/local/include/eo");
                                 
  @{$params{"make lib"}} = ( "-lm", "-lstdc++", "-L $HOME/usr/lib/",
                             "-llamarck", "-latat", "-ltinyxml",
                             "-lopt++", "-lnewmat",
                             "-lga", "-leoutils", "-leo" );
  $params{"CC"}  = "gcc";
  $params{"CXX"} = "g++";
  $params{"LD"}  = "g++";
  $params{"F77"}  = "g77";
  $params{"CXXFLAGS"}  = "-malign-double";
}
if ( $computer =~ /lester/ )
{
  @{$params{"make include"}} = ( "/usr/local/include", "$HOME/usr/include");
  @{$params{"make lib"}} = ( "-lm", "-lstdc++", "-L $HOME/usr/lib/",
                             "-llamarck", "-ltinyxml", "-lopt++", "-lnewmat" );
  $params{"CC"}  = "pgcc";
  $params{"CXX"} = "pgCC";
  $params{"LD"}  = "pgCC";
  $params{"F77"}  = "pgf90";
  $params{"CXXFLAGS"}  = "";

   
  $params{"cleanall"}[0]  = "- rm -f \$(addsuffix .ii, \$(basename \${OBJS})";
  $params{"cleanall"}[1]  = "- rm -f \$(addsuffix .ti, \$(basename \${OBJS})";

  push @{$params{"defs"}}, ( "__PORTLAND" ); 
}


# find source files in present directory
my %dependencies;
foreach $file ( <*.cc> )
{
  my $key = $file; $key =~ s/\.cc//;
  if ( $key =~ /\// )
  {
    $_ = $key; /(\S+)\/(\S+)/; $key = $2;
    $dependencies{$key}{"location"} = $1;
  }
  else
  {
    $dependencies{$key}{"location"} = ".";
  }
  if ( -e "$dependencies{$key}{'location'}/$key" . ".cc" )
    { $dependencies{$key}{"source"} = 1; }
  if ( -e "$dependencies{$key}{'location'}/$key" . ".h" )
  {
    $dependencies{$key}{"header"} = 1;
    push @{$dependencies{$key}{"depends on"}}, ( $key );
  }
}
# find orphaned header file in present directory
foreach $file ( <*.h> )
{
  if ( $file =~ /lada/ )
    { next; }
  if ( $file =~ /fftw_interface/ )
    { next; }
  my $key = $file; $key =~ s/\.h//;
  if ( ! (exists $dependencies{$key}{"header"}) )
  {
    $dependencies{$key}{"location"} = ".";
    $dependencies{$key}{"header"} = 1; 
  }
}

while ( get_dependencies() == 0 ) {};

write_make_file();
write_dependencies();


exit;

# now finds dependencies for each key
# if a new key is added, then its dependencies are also investigated
sub get_dependencies()
{
  my $key;
  my $new_stuff = 1;

  foreach $key ( keys %dependencies )
  {
    if ( exists $dependencies{$key}{"header"} ) 
    {
      open IN, "<$dependencies{$key}{'location'}/$key.h"
        or die "file $dependencies{$key}{'location'}/$key.h not found.\n";
    
      while ( <IN> )
      {
        if (/\#include(\s+|)(\"|\<)(\S+\/|)(\S+)\.h(\"|\>)/)
        {
          my $new_key = $4;
          if ( !( exists $dependencies{$new_key} ) )
          {
            foreach $location (  @{$params{"Includes"}} ) 
            { 
              if ( -e "$location/$new_key.h" )
              { 
                $dependencies{$new_key}{"location"} = $location;
                $dependencies{$new_key}{"header"} = 1;
                if ( -e "$location/$new_key.cc" )
                  { $dependencies{$new_key}{"source"} = 1; }
                last;
              }
            }
            if( exists $dependencies{$new_key} ) 
              { $new_stuff = 0; }

          } # end if ( !( exists $dependencies{$new_key} ) )
          
          add_to_dependencies( $key, $new_key ); 

        } # end if (/\#include(\s+|)\"(\S+)\.h\"/)
      } # end while ( <IN> )
      close IN;
    } # end if ( exists $dependencies{$key}{"header"} ) 
    if ( exists $dependencies{$key}{"source"} ) 
    {
      open IN, "<$dependencies{$key}{'location'}/$key.cc"
        or die "file $dependencies{$key}{'location'}/$key.cc not found.\n";
      while ( <IN> )
      {
        if (/\#include(\s+|)(\"|\<)(\S+\/|)(\S+)\.h(\"|\>)/)
        {
          my $new_key = $4;
          if ( !( exists $dependencies{$new_key} ) )
          {
            foreach $location (  @{$params{"Includes"}} ) 
            { 
              if ( -e "$location/$new_key.h" )
              { 
                $dependencies{$new_key}{"location"} = $location;
                $dependencies{$new_key}{"header"} = 1;
                if ( -e "$location/$new_key.cc" )
                  { $dependencies{$new_key}{"source"} = 1; }
                last;
              }
            }
            if( exists $dependencies{$new_key} ) 
              { $new_stuff = 0; }
           
          } # end if ( !( exists $dependencies{$new_key} ) )

          add_to_dependencies( $key, $new_key );

        } # end if (/\#include(\s+|)\"(\S+)\.h\"/)
      } # end while (<IN>)
      close IN;
    } # end  ( exists $dependencies{$key}{"source"} ) 
  }

  return $new_stuff;
}

sub template()
{
  open OUT, ">make_file_template" or die;
  printf OUT "CC     := $params{'CC'}\nCXX    := $params{'CXX'}\nF77    := $params{'F77'}\nLD     := $params{'LD'}\n";
  printf OUT "AR     := ar rc\nRANLIB := ranlib\nDEBUG = NO\n";
  printf OUT "\nFFLAGS     := \n";
  printf OUT  "CFLAGS     := \n";
  printf OUT  "CXXFLAGS   := $params{'CXXFLAGS'} \n";
  printf OUT "DEFS       := ";
  foreach $def ( @{$params{"defs"}} )
    { print OUT "-D ", $def, " "; }
  print OUT "\n\n";


  printf OUT "\nDEBUG_CFLAGS     := -Wall -Wno-format -g -O0 -Wno-unknown-pragmas -fbounds-check\n";
  printf OUT "DEBUG_CXXFLAGS    := \${DEBUG_CFLAGS} -D_DEBUG_LADA_ \n";
  printf OUT "DEBUG_LDFLAGS     := -g \n";
  printf OUT "\nRELEASE_CXXFLAGS := \${RELEASE_CFLAGS}\n";
  printf OUT "RELEASE_LDFLAGS   := \n";
  printf OUT "RELEASE_CFLAGS    := -Wall -Wno-unknown-pragmas -Wno-format -O3\n";
  printf OUT "\nLIBS:=  ";
  foreach $include (@{$params{"make lib"}})
    { print OUT " $include "; } 
  print OUT "\n";
  printf OUT "INCS := ";
  foreach $lib (@{$params{"make include"}})
    { print OUT "-I ", $lib, " "; } 
  print OUT "\n";
  printf OUT "\n\nifeq (YES, \${DEBUG})\n";
  printf OUT "   CFLAGS      := \${CFLAGS} \${DEBUG_CFLAGS} \n";
  printf OUT "   CXXFLAGS    := \${CXXFLAGS} \${DEBUG_CXXFLAGS}\n";
  printf OUT "   LDFLAGS     := \${LDFLAGS} \${DEBUG_LDFLAGS}\n";
  printf OUT "else\n";
  printf OUT "   CFLAGS      := \${CFLAGS} \${RELEASE_CFLAGS} -O3\n";
  printf OUT "   CXXFLAGS    := \${CXXFLAGS} \${RELEASE_CXXFLAGS} -O3 \n";
  printf OUT "   LDFLAGS     := \${LD_FLAGS} \${RELEASE_LDFLAGS} -O3\n";
  printf OUT "endif\n";
  printf OUT "\nCFLAGS   := \${CFLAGS}   \${DEFS}\n";
  printf OUT "CXXFLAGS := \${CXXFLAGS} \${DEFS}\n";
  printf OUT "\nOUTPUT := lada\n";
  printf OUT "\nall: \${OUTPUT} \n";
  printf OUT "\n\?SRCS :=\?\n";
  printf OUT "\nOBJS := \$(addsuffix .o,\$(basename \${SRCS}))\n";
  printf OUT "\n.PHONY: clean cleanall\n";
  printf OUT "\n\${OUTPUT}: \${OBJS} \n";
  printf OUT "\t\${LD} \${LDFLAGS} -o \$@ \${OBJS} \${LIBS} \${EXTRALIBS}\n";
  printf OUT "\n\${OBJS} : ${OBJSRCS}\n";
  printf OUT "\t\${CXX} -c \${CXXFLAGS} \${INCS} \$< -o \$@\n\n";
  printf OUT "\?dependencies\?";

            

  printf OUT "\n\nclean:\n\t- rm -f \${OBJS}\n";
  printf OUT "\t- rm -f \${OUTPUT}\n";
  if (exists $params{'cleanall'} )
  { 
    foreach $clean ( (@{$params{"cleanall"}}) )
    {
      printf OUT "\t%s\n", $clean;
    }
  }

  printf OUT "\n\ninclude .dependencies\n\n"; 

  close OUT;
}



sub  write_make_file()
{
  my @sorted_keys = sort{ sort_hash($a, $b) } keys %dependencies;
  template();

  open IN, "make_file_template" or die;
  open OUT, ">makefile" or die;
  
  while ( ($_=<IN> ) )
  {
    if ( /\?SRCS :=\?/ )
    {
      print OUT 'SRCS := ';
      my $i = 0;
      foreach $key ( @sorted_keys )
      {
        if ( exists $dependencies{$key}{"source"} )
        {
          if ( $i % 5 == 0 and $i != 0 )
            { print OUT "\\\n\t"; }

          if ( $dependencies{$key}{"location"} eq "." )
            { print OUT $key, ".cc "; $i++; }
          else
            { print OUT $dependencies{$key}{"location"},"/", $key, ".cc "; $i++; }
        }
      }
      print OUT "\n";
      next;
    }
    elsif (/\?dependencies\?/)
    { next; }
    print OUT $_;
  }
  
  close IN;
  close OUT;
  system "rm make_file_template";

}

sub add_to_dependencies(\$\$)
{
  my $key = $_[0];
  my $new_key = $_[1];

  if ( ! exists $dependencies{$key}{"depends on"} )
  {
    push @{$dependencies{$key}{"depends on"}}, ( "$new_key" ); 
    return;
  }

  $w = 0;
  foreach $dep ( @{$dependencies{$key}{"depends on"}} )
  {
     if ( $dep eq $new_key )
       { last; }
     $w ++ ;
  }
  
  if ( $w == scalar( @{$dependencies{$key}{"depends on"}} ) )
    { push @{$dependencies{$key}{"depends on"}}, ( "$new_key" ); }
}


sub write_dependencies()
{
  open OUT, ">.dependencies"
    or die " could not open .dependencies for writing \n";

  my @sorted_keys = sort{ sort_hash($a, $b) } keys %dependencies;
  foreach $key ( @sorted_keys )
  {    
    if ( not exists $dependencies{$key}{"source"} )
      { next; }
    if ( $dependencies{$key}{"location"} ne "." )
      { print OUT $dependencies{$key}{"location"}, "/", $key, ".o: ",
                  $dependencies{$key}{"location"}, "/", $key, ".cc ";}
    else
      { print OUT $key, ".o: ", $key, ".cc ";}
    my $i = 0;
    foreach $dep ( sort { sort_hash($a,$b) } @{$dependencies{$key}{"depends on"}} )
    {
      if ( $dep =~ $key )
        { next; }
      if( $i % 1 == 0 and $i != 0)
        { print OUT " \\\n\t"; }
      if ( (exists $dependencies{$dep}{"source"}) and
           (exists $dependencies{$dep}{"header"})     )
      {
        if ( $dependencies{$dep}{"location"} ne "." )
          { print OUT $dependencies{$dep}{"location"}, "/"; }
        print OUT $dep, ".o "; $i++;
        next;
      }
      if ( (exists $dependencies{$dep}{"source"}) )
      {
        if ( $dependencies{$dep}{"location"} ne "." )
          { print OUT $dependencies{$dep}{"location"}, "/"; }
        print OUT $dep, ".cc "; $i++;
      }
      if ( exists $dependencies{$dep}{"header"} )
      {
        if ( exists $dependencies{$dep}{"source"} )
          { print OUT " \\\n\t"; }
        if ( $dependencies{$dep}{"location"} ne "." )
          { print OUT $dependencies{$dep}{"location"}, "/"; }
        print OUT $dep, ".h "; $i++;
      }
      
    }
    print OUT "\n\n";
  }    

}

sub sort_hash ()
{
  my $a = $_[0], $b = $_[1];

  if(     ($dependencies{$a}{"location"} eq ".")
      and ($dependencies{$b}{"location"} ne ".") )
    { return 0; }
  if(     ($dependencies{$a}{"location"} ne ".")
      and ($dependencies{$b}{"location"} eq ".") )
    { return 1; }
  if( $dependencies{$a}{"location"} ne $dependencies{$b}{"location"} )
    { return ($dependencies{$a}{"location"} > $dependencies{$b}{"location"} ); } 
  if ($a > $b)
    { return 0; }

  return 1;
}

sub is_in_compile_location($)
{
  my $location = $_[0];
  my $result = 0;
  foreach $cloc ( @{$params{"compile"}} )
  {
    if ( $location eq $cloc )
    {
      $result = 1 ;
      last;
    }
  }
  return $result;
}

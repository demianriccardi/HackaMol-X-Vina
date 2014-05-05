HackaMol-X-Vina
===============
HackaMol extension for running Autodock Vina  

VERSION
========
developer version 0.00_4 
Available for testing from cpan.org:

please see *[HackaMol::X::Vina on MetaCPAN](https://metacpan.org/release/DEMIAN/HackaMol-X-Vina-0.00_4) for formatted documentation.

SYNOPSIS
============
         use Modern::Perl;
         use HackaMol;
         use HackaMol::X::Vina;
         use Math::Vector::Real;
         
         my $receptor = "receptor.pdbqt";
         my $ligand   = "lig.pdbqt",
         my $rmol     = HackaMol -> new( hush_read=>1 ) -> read_file_mol( $receptor ); 
         my $lmol     = HackaMol -> new( hush_read=>1 ) -> read_file_mol( $ligand ); 
         my $fh = $lmol->print_pdb("lig_out.pdb");
  
         my @centers = map  {$_ -> xyz}
                       grep {$_ -> name    eq "OH" }
                       grep {$_ -> resname eq "TYR"} $rmol -> all_atoms;
         
         foreach my $center ( @centers ){
         
             my $vina = HackaMol::X::Vina -> new(
                 receptor       => $receptor,
                 ligand         => $ligand,
                 center         => $center,
                 size           => V( 20, 20, 20 ),
                 cpu            => 4,
                 exhaustiveness => 12,
                 exe            => '~/bin/vina',
                 scratch        => 'tmp',
             );
         
             my $mol = $vina->dock_mol(3); # fill mol with 3 binding configurations 
         
             printf ("Score: %6.1f\n", $mol->get_score($_) ) foreach (0 .. $mol->tmax);          

             $mol->print_pdb_ts([0 .. $mol->tmax], $fh); 

           }

           $_->segid("hgca") foreach $rmol->all_atoms; #for vmd rendering cartoons.. etc
           $rmol->print_pdb("receptor.pdb");

DESCRIPTION
============
HackaMol::X::Vina provides an interface to AutoDock Vina. This class does not include the AutoDock Vina program, which is 
[released under a very permissive Apache license](http://vina.scripps.edu/manual.html#license), with few restrictions on 
commercial or non-commercial use, or on the derivative works, such is this. Follow these 
[instructions ] (http://vina.scripps.edu/manual.html#installation) to acquire the program. Most importantly, if you use this 
interface effectively, please be sure to cite AutoDock Vina in your work:

O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461 

Since HackaMol has no pdbqt writing capabilities (yet, HackaMol can read pdbqt files), the user is required to provide 
those  files. This is still a work in progress and the API may still change. Documentation will improve as API 
gets more stable... comments welcome! 


#!/usr/bin/env perl

use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Fatal qw(lives_ok dies_ok);
use Test::Dir;
use Test::Warn;
use HackaMol::X::Vina;
use HackaMol;
use Math::Vector::Real;
use File::chdir;
use Path::Tiny;
use Cwd;

my $cwd      = getcwd;
my $receptor = path('t/lib/receptor.pdbqt'); 
my $ligand1  = path('t/lib/ligand1.pdbqt' ); 
my $ligand2  = path('t/lib/ligand2.pdbqt' ); 
my $rmol     = HackaMol->new(hush_read=>1)->read_file_mol($receptor);    
my $lmol1    = HackaMol->new(hush_read=>1)->read_file_mol($ligand1);    
my $lmol2    = HackaMol->new(hush_read=>1)->read_file_mol($ligand2);    

my $obj = HackaMol::X::Vina->new(
        receptor       => $receptor->basename,
        ligand         => $ligand1->basename,
        in_fn          => "conf.txt",
        center         => V(0,1,2),
        size           => V(20,20,20),
        cpu            => 1,
        num_modes      => 2,
        exhaustiveness => 8,
        exe            => '~/bin/vina',
        scratch        => 't/tmp',
        seed           => 314159,
    );

my $input = $obj->map_input;
print $input;

done_testing();


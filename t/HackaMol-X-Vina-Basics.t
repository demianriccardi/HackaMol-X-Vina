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
use Cwd;

BEGIN {
    use_ok('HackaMol::X::Vina');
}

my $cwd = getcwd;

# coderef

{    # test HackaMol class attributes and methods

    my @attributes = qw(mol map_in map_out);
    my @methods    = qw(map_input map_output);

    my @roles = qw(HackaMol::ExeRole HackaMol::PathRole);

    map has_attribute_ok( 'HackaMol::X::Vina', $_ ), @attributes;
    map can_ok( 'HackaMol::X::Vina', $_ ), @methods;
    map does_ok( 'HackaMol::X::Vina', $_ ), @roles;

}

my $mol = HackaMol::Molecule->new();
my $obj;

{    # test basic functionality

    lives_ok {
        $obj = HackaMol::X::Vina->new(
            receptor => 't/lib/receptor.pdbqt',
            ligand   => 't/lib/receptor.pdbqt',
            center => V( 0,  1,  2 ),
            size   => V( 10, 11, 12 ),
        );
    }
    'creation without required mol lives';

    is( $obj->center_x, 0,  "center_x" );
    is( $obj->center_y, 1,  "center_y" );
    is( $obj->center_z, 2,  "center_z" );
    is( $obj->size_x,   10, "size_x" );
    is( $obj->size_y,   11, "size_y" );
    is( $obj->size_z,   12, "size_z" );

    lives_ok {
        $obj = HackaMol::X::Vina->new( 
                                      mol => $mol , 
                                      receptor => 't/lib/receptor.pdbqt', 
                                      ligand   => 't/lib/ligand.pdbqt', 
                                     );
    }
    'creation of an obj with mol';

    dir_not_exists_ok( "t/tmp", 'scratch directory does not exist yet' );

    lives_ok {
        $obj = HackaMol::X::Vina->new( 
                                      mol => $mol, 
                                      exe => "vina",
                                      receptor => 't/lib/receptor.pdbqt', 
                                      ligand   => 't/lib/ligand.pdbqt', 
                                     );
    }
    'creation of an obj with exe';

    dir_not_exists_ok( "t/tmp", 'scratch directory does not exist yet' );

    is($obj->in_fn, 'conf.txt', "default configuration file conf.txt" );

    is(
        $obj->command,
        $obj->exe . " --config " . $obj->in_fn,
        "command set to exe and input"
    );

    lives_ok {
        $obj = HackaMol::X::Vina->new(
            mol     => $mol,
            exe     => "vina",
            in_fn   => "foo.inp",
            receptor => 't/lib/receptor.pdbqt', 
            ligand   => 't/lib/ligand.pdbqt', 
            scratch => "t/tmp"
        );
    }
    'Test creation of an obj with exe in_fn and scratch';

    dir_exists_ok( $obj->scratch, 'scratch directory exists' );
    is(
        $obj->command,
        $obj->exe . " --config " . $obj->in_fn,
        "command set to exe and input"
    );
    is( $obj->scratch, "$cwd/t/tmp", "scratch directory" );

    lives_ok {
        $obj = HackaMol::X::Vina->new(
            mol     => $mol,
            exe     => "vina",
            in_fn   => "foo.inp",
            scratch => "t/tmp",  
            receptor => 't/lib/receptor.pdbqt', 
            ligand   => 't/lib/ligand.pdbqt', 
            command => "nonsense",
        );
    }
    'test building of an obj with exisiting scratch  and command attr';

    is( $obj->command, "nonsense",
        "command attr not overwritten during build" );

    $obj->command( $obj->build_command );
    is( $obj->command, $obj->exe . " --config " . $obj->in_fn,
        "command reset" );

    $obj->scratch->remove_tree;
    dir_not_exists_ok( "t/tmp", 'scratch directory deleted' );

    lives_ok {
        $obj = HackaMol::X::Vina->new(
            mol        => $mol,
            exe        => "vina",
            in_fn      => "foo.inp",
            scratch    => "t/tmp",
            out_fn     => "foo.out",
            command    => "nonsense",
            exe_endops => "tackon",
            receptor => 't/lib/receptor.pdbqt', 
            ligand   => 't/lib/ligand.pdbqt', 
        );
    }
    'test building of an obj with out_fn';

    $obj->command( $obj->build_command );
    is(
        $obj->command,
        $obj->exe . " --config " . $obj->in_fn->stringify,
        "big command ignores redirect to output"
    );

    $obj->scratch->remove_tree;
    dir_not_exists_ok( "t/tmp", 'scratch directory deleted' );

}

{    # test the map_in and map_out

    $obj = HackaMol::X::Vina->new(
        mol            => $mol,
        receptor       => 't/lib/receptor.pdbqt', 
        ligand         => 't/lib/ligand.pdbqt', 
        in_fn          => "foo.inp",
        center         => V( 0, 1, 2 ),
        size           => V( 20, 20, 20 ),
        cpu            => 4,
        num_modes      => 1,
        exhaustiveness => 12,
        exe            => '~/bin/vina',
        scratch        => 't/tmp',
        homedir        => '.',
    );

    my $input = $obj->map_input; 
    $CWD = $obj->scratch;
    my $input2 = $obj->in_fn->slurp;
    is( $input, $input2,
        "input written to scratch is that returned by map_input" );
    $CWD = $obj->homedir;
    $obj->scratch->remove_tree;
    dir_not_exists_ok( "t/tmp", 'scratch directory deleted' );

}

done_testing();


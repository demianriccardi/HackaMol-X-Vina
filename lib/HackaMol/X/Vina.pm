  package HackaMol::X::Vina;
  #ABSTRACT: HackaMol extension for running Autodock Vina  
  use Moose;
  use Moose::Util::TypeConstraints;
  use Carp;
  use namespace::autoclean;
 
  with qw(HackaMol::ExeRole HackaMol::PathRole);

  has 'mol' => (
    is  => 'ro',
    isa => 'HackaMol::Molecule',
  );

  has 'receptor'  => (is => 'ro', isa => 'Str', predicate => 'has_receptor');
  has 'ligand'    => (is => 'ro', isa => 'Str', predicate => 'has_ligand');
  has 'cpu'       => (is => 'ro', isa => 'Int', predicate => 'has_cpu');
  has 'num_modes' => (is => 'ro', isa => 'Int', predicate => 'has_num_modes');

  has $_ => (
      is => 'ro', isa => 'Num', predicate => "has_$_",
  ) foreach qw(center_x center_y center_z);

  has $_ => (
    is => 'ro', isa => 'Int', predicate => "has_$_",
  ) foreach qw(energy_range exhaustiveness);

  sub BUILD {
    my $self = shift;

    if ( $self->has_scratch ) {
        $self->scratch->mkpath unless ( $self->scratch->exists );
    }
    $self->exe("~/bin/vina") unless $self->has_exe;

    return;
  }

  sub build_command {
    my $self = shift;
    my $cmd;
    $cmd  = $self->exe;
    $cmd .= " " . $self->in_fn->stringify    if $self->has_in_fn;
    $cmd .= " " . $self->exe_endops          if $self->has_exe_endops;
    $cmd .= " > " . $self->out_fn->stringify if $self->has_out_fn;

    # no cat of out_fn because of options to run without writing, etc
    return $cmd;
}

  sub write_input {
    my $self  = shift;
    local $CWD = $self->scratch if ( $self->has_scratch );

    my $input ;
    $input   .= "out       = "    . $self->out_fn->stringify . "\n" if $self->has_out;
    $input   .= "cpu       = "    . $self->cpu               . "\n" if $self->has_cpu;
    $input   .= "num_modes = "    . $self->num_modes         . "\n" if $self->has_num_modes;
    $input   .= "log       = "    . $self->err->stringify    . "\n" if $self->has_err;
    foreach my $cond (qw(receptor ligand cpu num_modes energy_range exhaustiveness)) {
      my $condition = "has_$cond";
      $input .= "$cond = ". $self->$cond . "\n" if $self->$condition;
    }
    foreach my $metric (qw(center_x center_y center_z size_x size_y size_z)) {
      $input .= "$metric = ". $self->$metric . "\n";
    }
  
    $self->in_fn->spew($input);
      
  }

  __PACKAGE__->meta->make_immutable;

  1;

__END__

=head1 SYNOPSIS

  


  package HackaMol::X::Vina;
  #ABSTRACT: HackaMol extension for running Autodock Vina  
  use Moose;
  use MooseX::StrictConstructor;
  use Moose::Util::TypeConstraints;
  use Math::Vector::Real;
  use namespace::autoclean;
  use Carp;

  with qw(HackaMol::X::ExtensionRole);

  has 'receptor'  => (is => 'ro', isa => 'Str', predicate => 'has_receptor');
  has 'ligand'    => (is => 'ro', isa => 'Str', predicate => 'has_ligand');
  has 'cpu'       => (is => 'ro', isa => 'Int', predicate => 'has_cpu');
  has 'num_modes' => (is => 'ro', isa => 'Int', predicate => 'has_num_modes');

  has $_ => (
      is => 'rw', isa => 'Num', predicate => "has_$_",
  ) foreach qw(center_x center_y center_z size_x size_y size_z);

  has $_ => (
      is => 'ro', isa => 'Math::Vector::Real', predicate => "has_$_",
  ) foreach qw(center size);

  has $_ => (
    is => 'ro', isa => 'Int', predicate => "has_$_",
  ) foreach qw(energy_range exhaustiveness seed);


  sub BUILD {
    my $self = shift;

    if ( $self->has_scratch ) {
        $self->scratch->mkpath unless ( $self->scratch->exists );
    }

#  setting center via MVR overrides whatever is set via center_x etc.
    if ($self->has_center){
      my ($x,$y,$z) = @{$self->center};
      $self->center_x($x);  
      $self->center_y($y);  
      $self->center_z($z);  
    }

    if ($self->has_size){
      my ($x,$y,$z) = @{$self->size};
      $self->size_x($x); 
      $self->size_y($y);
      $self->size_z($z);
    }
    
    unless ( $self->has_command ) {
        return unless ( $self->has_exe );
        my $cmd = $self->build_command;
        $self->command($cmd);
    }
    return;
  }

  #required methods
  sub build_command {
    my $self = shift;
    my $cmd;
    $cmd  = $self->exe;
    $cmd .= " --config " . $self->in_fn->stringify  if $self->has_in_fn;
    # we always capture output 
    return $cmd;
  }

  sub _build_map_in{
    
    my $sub_cr = sub { 
                      my $self = shift;
                      $self->write_input;
                     };
    return $sub_cr;
  }

  sub _build_map_out{
    my $sub_cr = sub { return (@_) };
    return $sub_cr;
  }

  sub write_input {
    my $self  = shift;

    unless ($self->has_in_fn) {
      croak "no vina in_fn for writing input";
    }

    my $input ;
    $input   .= "out       = "    . $self->out_fn->stringify . "\n"    if $self->has_out_fn;
    $input   .= "cpu       = "    . $self->cpu               . "\n"    if $self->has_cpu;
    $input   .= "num_modes = "    . $self->num_modes         . "\n"    if $self->has_num_modes;
    $input   .= "log       = "    . $self->log_fn->stringify    . "\n" if $self->has_log_fn;
    foreach my $cond (qw(receptor ligand cpu num_modes energy_range exhaustiveness seed)) {
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

  


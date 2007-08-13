#!/usr/bin/perl -w
#
#
###############################################################

package WriteXML;

=pod

=head1 NAME

WriteXML

=head1 AUTHOR

Kelvin Li

=head1 DESCRIPTION

Function calls to generate XML

=cut

use Carp;
use strict;

###############################################################################

sub element{

=head2 str element(str name, hash attributes, int indent, str data);

This subroutine will generate an XML element.

<name> is the name of the element

<attributes> is a hash that contains the key/value of the attributes

<indent> is a boolean that tells whether to indent the contents of the element

<data> is the data that should go inside of the element.

If you are not sure if your data has & or <, then you should apply the 
escape function to your string.

=cut

	my $name=shift;
	my $attributes=shift;
	my $indent=shift;
	my $data=shift;

	my $xml;

	# Start the element
	$xml="<$name";

	# Output the attributes
	foreach my $key(sort keys %{$attributes}){
		$xml .= " $key=\"${$attributes}{$key}\"";
	}
	
	if(!defined($data) || $data eq ""){
		# If data is null, output the shortened element type
		$xml.="/>\n"
	}else{
		if($indent){
			# If using the indent option, indent all new lines
			my $last_data_char=chop $data;
			if($last_data_char ne "\n"){
				$data.=$last_data_char;
			}
			$data=~s/\n/\n\t/g;
			$data="\n\t".$data;
			$xml.=">$data\n</$name>\n";
		}else{
			#If not using the indent option, just slap the end tag on
			$xml.=">$data</$name>\n";
		}
	}

	return $xml;
}

###############################################################################

sub version{

=head2 str version(str version, str encoding);

Returns a XML version string with the version/encoding information specified.

=cut

	my $version=shift;
	my $encoding=shift;

	return "<?xml version=\"$version\" encoding=\"$encoding\"?>\n";

}

###############################################################################

sub comment{

=head2 str comment(str comment);

Returns a XML comment string with your comment inside

=cut

	my $comment=shift;

	return "<!--$comment-->\n";

}

###############################################################################

sub escape{

=head2 str escape(str data);

Returns the string you entered with the & and < escaped
with the &amp and &lt, respectively.

=cut

	my $data=shift;

	$data=~s/&/&amp/g;
	$data=~s/</&lt/g;
	
	return($data);

}

###############################################################################

sub wrapLongStr{

=head2 str wrapLongStr(str data, int width);

Returns the input string with line returns.
Each line becomes limited to width specified.
Useful for printing out sequences.

=cut
	my $data=shift;
	my $width=shift;
	my $out_data="";

	my $length=length($data);
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		$out_data.=substr($data, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
	return $out_data;

}

1;


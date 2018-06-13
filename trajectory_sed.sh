#!/bin/sh

sed -n '/         x   /		{
	N
	/\n[^*]\{35\}/			{
	N
	N
	N
	N
	N
//p
}}' $1 > $2

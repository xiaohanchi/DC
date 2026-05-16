#!/bin/bash

bjobs -u all -r | grep "short" | grep "xchi" | wc -l
bjobs -u all -r | grep "medium" | grep "xchi" | wc -l
bjobs -u all -r | grep "long" | grep "xchi" | wc -l

bjobs -u all -p | grep "short" | grep "xchi" | wc -l
bjobs -u all -p | grep "medium" | grep "xchi" | wc -l
bjobs -u all -p | grep "long" | grep "xchi" | wc -l

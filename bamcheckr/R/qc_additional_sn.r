###############################################################################
# qc_additional_sn.r
###############################################################################
#
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Sam Nicholls <sn8@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

qc_additional_sn <- function(bamcheck) {

    summary_numbers <- bamcheck$data$SN
    indel_distribution <- bamcheck$data$ID
    read_lengths <- bamcheck$data$RL
    insert_sizes <- bamcheck$data$IS

    duplicated_reads_percentage <- duplicated_reads_pct(as.numeric(summary_numbers[summary_numbers[1] == "reads mapped:"][2]),
                                                        as.numeric(summary_numbers[summary_numbers[1] == "reads duplicated:"][2]))
    insertions_to_deletions_ratio <- ins_to_del_ratio(indel_distribution)
    overlapping_base_duplicate_percentage <- overlapping_base_duplicate_pct(read_lengths, insert_sizes)

    outdata <- data.frame(
                    variable = c("duplicate.read.percentage:",
                                 "ins.to.del.ratio:",
                                 "overlapping.base.duplicate.percent:"
                    ),
                    value = c(duplicated_reads_percentage,
                              insertions_to_deletions_ratio,
                              overlapping_base_duplicate_percentage
                    ))
    return(outdata)
}

###############################################################################
duplicated_reads_pct <- function(mapped_reads, duplicated_reads){
    return (100.00 * duplicated_reads / mapped_reads)
}
#bases_mapped_pct  = 100 * $bases_mapped_c / $clip_bases;
#    return (100.00 * bases_mapped / clipped_bases)
#}

ins_to_del_ratio <- function(indel_distribution){
    return (sum(indel_distribution[2]) / sum(indel_distribution[3]))
}

overlapping_base_duplicate_pct <- function(read_lengths, insert_sizes){
    # If consistent read length
    if(nrow(read_lengths) == 1 && nrow(insert_sizes) > 0){

        seq_length <- read_lengths[1][1]

        insert_sizes$short_paired_reads <- apply(insert_sizes, 1, function(row, seq_length){
            is <- row[1]
            pairs_total <- row[2]
            if((seq_length * 2) > is){
                return(pairs_total)
            }
            else { return (0.0) }
        }, seq_length)

        insert_sizes$normal_paired_reads <- apply(insert_sizes, 1, function(row, seq_length){
            is <- row[1]
            pairs_total <- row[2]
            if(!(seq_length * 2) > is){
                return(pairs_total)
            }
            else { return (0.0) }
        }, seq_length)

        insert_sizes$total_paired_reads <- apply(insert_sizes, 1, function(row, seq_length){
            pairs_total <- row[2]
            return(pairs_total)
        }, seq_length)

        insert_sizes$dup_mapped_bases <- apply(insert_sizes, 1, function(row, seq_length){
            is <- row[1]
            pairs_total <- row[2]
            if((seq_length * 2) > is){
                return(as.numeric(pairs_total * ((seq_length * 2) - is)))
            }
            else { return (0.0) }
        }, seq_length)

        insert_sizes$tot_mapped_bases <- apply(insert_sizes, 1, function(row, seq_length){
            is <- row[1]
            pairs_total <- row[2]
            if((seq_length * 2) > is){
                return(as.numeric(pairs_total * (seq_length * 2)))
            }
            else { return (0.0) }
        }, seq_length)

        t_dup_mapped_bases <- sum(insert_sizes$dup_mapped_bases)
        t_tot_mapped_bases <- sum(insert_sizes$tot_mapped_bases)
        return (100.00 * t_dup_mapped_bases / t_tot_mapped_bases)
    }
    else{
        print("[WARN] Variable read lengths detected in overlapping_base_duplicate_pct")
        return
    }
}


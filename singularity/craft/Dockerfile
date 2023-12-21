FROM annadalmolin/craft:v1.0

COPY script_sequence_extraction.sh /scripts/
COPY pipeline_predictions.sh /scripts/
ENTRYPOINT []
CMD []

#run these commands to open up permissions within the tool dir and home dir
RUN chmod -vR 777 /tools
RUN chmod -vR 777 /home

#run this command to change paths for /data to /annadalmolin_craft_root
COPY files_for_edit.txt /root/
RUN cd /root && cat files_for_edit.txt | xargs -tI {} sed -r -i  "s/\/data/\/annadalmolin_craft_root/g" {} 


#run these commands to allow usage of the symlink
RUN rmdir /data/
RUN mkdir -v /containing_dir
RUN chmod -v 777 /containing_dir
RUN ln -vs /containing_dir/link_to_target /annadalmolin_craft_root


# RUN THIS TO COMPLETE THE LINK ONCE INSIDE THE CONTAINER
#ln -vs /mnt/path/to/target /containg_dir/link_to_target
# e.g. ln -vs ${PWD} /containing_dir/link_to_target

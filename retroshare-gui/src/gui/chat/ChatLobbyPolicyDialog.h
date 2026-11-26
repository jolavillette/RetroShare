/*******************************************************************************
 * retroshare-gui/src/gui/chat/ChatLobbyPolicyDialog.h                         *
 *                                                                             *
 * Copyright (C) 2026, Retroshare Team <retroshare.project@gmail.com>          *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Affero General Public License as              *
 * published by the Free Software Foundation, either version 3 of the          *
 * License, or (at your option) any later version.                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful,             *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
 * GNU Affero General Public License for more details.                         *
 *                                                                             *
 * You should have received a copy of the GNU Affero General Public License    *
 * along with this program. If not, see <https://www.gnu.org/licenses/>.       *
 *                                                                             *
 *******************************************************************************/
#ifndef CHATLOBBYPOLICYDIALOG_H
#define CHATLOBBYPOLICYDIALOG_H

#include <QDialog>

#include <retroshare/rschats.h>

class QCheckBox;
class QTreeWidget;

/// A hosted room: a chat room whose relay rules a node, the host, holds. On
/// the host, or on a node
/// that was told of no host: the two switches and the participant table,
/// where the admins are chosen. On a participant of a room another node
/// hosts: what the host told it, read-only.
class ChatLobbyPolicyDialog : public QDialog
{
	Q_OBJECT

public:
	ChatLobbyPolicyDialog(const ChatLobbyId& lobbyId, QWidget *parent = nullptr);

private slots:
	void apply();
	void enableChanged(bool enabled);

private:
	void buildEditor();
	void buildRoleView(const ChatLobbyRole& role);
	void fillParticipants(const ChatLobbyPolicy& policy);

	ChatLobbyId mLobbyId;
	ChatLobbyPolicy mPolicy;

	QCheckBox *mEnabled;
	QCheckBox *mMembersSee;
	QCheckBox *mMembersTalk;
	QTreeWidget *mParticipants;
};

#endif
